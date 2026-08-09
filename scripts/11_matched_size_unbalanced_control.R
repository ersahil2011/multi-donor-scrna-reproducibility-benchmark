# ============================================================
# 11_matched_size_unbalanced_control.R
#
# Purpose:
# Compare donor-balanced sampling with a matched-size control
# that preserves the original unequal donor proportions.
#
# Original RDS objects will only be read, never overwritten.
# ============================================================

rm(list = ls())
gc()

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "mclust",
  "aricode",
  "dplyr",
  "tidyr",
  "ggplot2",
  "future"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0) {
  stop(
    "Missing packages: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(mclust)
  library(aricode)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(future)
})

future::plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

cat("Step 2 packages loaded successfully.\n")


# ============================================================
# 2. Define files and load source objects
# ============================================================

output_dir <- file.path(
  "reviewer_revision_outputs",
  "11_matched_size_unbalanced_control"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

input_files <- c(
  post_filter_object = "sce_after_doublet_removal.rds",
  reference_object =
    "sce_final_clusters_locked_rpca_pcs20_res06.rds"
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {
  stop(
    "Missing input files: ",
    paste(missing_files, collapse = ", ")
  )
}

cat("Loading post-doublet-removal object...\n")
post_obj <- readRDS(
  unname(input_files["post_filter_object"])
)

cat("Loading full-data reference object...\n")
reference_obj <- readRDS(
  unname(input_files["reference_object"])
)

stopifnot(
  inherits(post_obj, "Seurat"),
  inherits(reference_obj, "Seurat"),
  identical(dim(post_obj), dim(reference_obj)),
  identical(colnames(post_obj), colnames(reference_obj))
)

DefaultAssay(post_obj) <- "RNA"

cat(
  "Objects loaded and aligned:",
  nrow(post_obj), "genes and",
  ncol(post_obj), "cells.\n"
)

cat("Output folder:", output_dir, "\n")



# ============================================================
# 3. Define the matched-size unbalanced sampling design
# ============================================================

donor_column <- "sample_id"

if (!donor_column %in% colnames(post_obj@meta.data)) {
  stop("The sample_id donor column is missing.")
}

donor_counts <- table(
  post_obj@meta.data[[donor_column]]
)

if (length(donor_counts) != 5) {
  stop("Exactly five donors were expected.")
}

cells_per_donor_balanced <- min(
  as.integer(donor_counts)
)

matched_total_cells <-
  cells_per_donor_balanced * length(donor_counts)

sampling_design <- data.frame(
  donor = names(donor_counts),
  available_cells = as.integer(donor_counts),
  stringsAsFactors = FALSE
)

sampling_design$original_proportion <-
  sampling_design$available_cells /
  sum(sampling_design$available_cells)

sampling_design$unrounded_target <-
  sampling_design$original_proportion *
  matched_total_cells

sampling_design$target_cells <-
  floor(sampling_design$unrounded_target)

remaining_cells <-
  matched_total_cells -
  sum(sampling_design$target_cells)

if (remaining_cells > 0) {
  fractional_remainder <-
    sampling_design$unrounded_target -
    sampling_design$target_cells
  
  receive_extra_cell <- order(
    fractional_remainder,
    decreasing = TRUE
  )[seq_len(remaining_cells)]
  
  sampling_design$target_cells[
    receive_extra_cell
  ] <- sampling_design$target_cells[
    receive_extra_cell
  ] + 1L
}

sampling_design$matched_sample_proportion <-
  sampling_design$target_cells /
  sum(sampling_design$target_cells)

stopifnot(
  sum(sampling_design$target_cells) ==
    matched_total_cells,
  all(
    sampling_design$target_cells <=
      sampling_design$available_cells
  ),
  all(sampling_design$target_cells > 0)
)

write.csv(
  sampling_design,
  file.path(
    output_dir,
    "Table_matched_size_unbalanced_sampling_design.csv"
  ),
  row.names = FALSE
)

print(sampling_design, row.names = FALSE)

cat(
  "Matched-size control total:",
  matched_total_cells,
  "cells.\n"
)




# ============================================================
# 4. Identify donor-specific raw count layers
# ============================================================

rna_layer_names <- Layers(
  post_obj[["RNA"]]
)

count_layer_names <- rna_layer_names[
  grepl(
    "^counts($|\\.)",
    rna_layer_names,
    ignore.case = TRUE
  )
]

print(count_layer_names)

if (length(count_layer_names) != 5) {
  stop(
    "Expected five raw count layers, but detected ",
    length(count_layer_names), "."
  )
}

cat("Five raw count layers verified.\n")






# ============================================================
# 5. Combine and align raw count layers
# ============================================================

feature_universe <- rownames(post_obj)

extract_aligned_count_layer <- function(layer_name) {
  
  cat("Extracting:", layer_name, "\n")
  
  layer_matrix <- LayerData(
    object = post_obj,
    assay = "RNA",
    layer = layer_name
  )
  
  unexpected_features <- setdiff(
    rownames(layer_matrix),
    feature_universe
  )
  
  if (length(unexpected_features) > 0) {
    stop(
      "Unexpected features detected in layer: ",
      layer_name
    )
  }
  
  missing_features <- setdiff(
    feature_universe,
    rownames(layer_matrix)
  )
  
  if (length(missing_features) > 0) {
    
    zero_matrix <- Matrix::Matrix(
      0,
      nrow = length(missing_features),
      ncol = ncol(layer_matrix),
      sparse = TRUE,
      dimnames = list(
        missing_features,
        colnames(layer_matrix)
      )
    )
    
    layer_matrix <- rbind(
      layer_matrix,
      zero_matrix
    )
  }
  
  layer_matrix <- layer_matrix[
    feature_universe,
    ,
    drop = FALSE
  ]
  
  cat(
    "Aligned dimensions:",
    paste(dim(layer_matrix), collapse = " x "),
    "\n"
  )
  
  layer_matrix
}

count_matrix_list <- lapply(
  count_layer_names,
  extract_aligned_count_layer
)

raw_counts <- do.call(
  cbind,
  count_matrix_list
)

if (anyDuplicated(colnames(raw_counts)) > 0) {
  stop("Duplicate cell barcodes detected.")
}

common_cells <- colnames(post_obj)[
  colnames(post_obj) %in% colnames(raw_counts)
]

if (length(common_cells) != ncol(post_obj)) {
  stop("Not all post-filtered cells were recovered.")
}

raw_counts <- raw_counts[
  ,
  common_cells,
  drop = FALSE
]

metadata_source <- post_obj@meta.data[
  common_cells,
  ,
  drop = FALSE
]

aligned_donor_counts <- table(
  metadata_source[[donor_column]]
)

stopifnot(
  identical(
    as.integer(aligned_donor_counts),
    as.integer(donor_counts)
  )
)

cat(
  "Combined raw count dimensions:",
  paste(dim(raw_counts), collapse = " x "),
  "\n"
)

print(aligned_donor_counts)



# ============================================================
# 6. Define reference labels and control iterations
# ============================================================

reference_cluster_column <- "cluster_final"

reference_cluster_data <- FetchData(
  object = reference_obj,
  vars = reference_cluster_column,
  clean = "none"
)

full_reference_vector <- as.character(
  reference_cluster_data[, 1]
)

names(full_reference_vector) <- rownames(
  reference_cluster_data
)

if (anyNA(full_reference_vector[common_cells])) {
  stop("Reference labels are missing for some aligned cells.")
}

donor_vector <- as.character(
  metadata_source[, donor_column]
)

names(donor_vector) <- rownames(metadata_source)

cells_by_donor <- split(
  common_cells,
  donor_vector[common_cells]
)

target_cells_by_donor <- setNames(
  sampling_design$target_cells,
  sampling_design$donor
)

target_cells_by_donor <- target_cells_by_donor[
  names(cells_by_donor)
]

stopifnot(
  !anyNA(target_cells_by_donor),
  sum(target_cells_by_donor) == matched_total_cells,
  all(
    target_cells_by_donor <=
      lengths(cells_by_donor)
  )
)

number_of_iterations <- 50
iteration_seeds <- 2001:2050

write.csv(
  data.frame(
    iteration = seq_len(number_of_iterations),
    seed = iteration_seeds
  ),
  file.path(
    output_dir,
    "matched_size_unbalanced_iteration_seeds.csv"
  ),
  row.names = FALSE
)

cat("Reference clustering verified.\n")
cat("Planned control iterations:", number_of_iterations, "\n")
cat("Control cells per iteration:", matched_total_cells, "\n")
print(target_cells_by_donor)






# ============================================================
# 7. Define one matched-size unbalanced RPCA iteration
# ============================================================

run_unbalanced_control_iteration <- function(
    iteration_number,
    iteration_seed
) {
  
  iteration_start <- unname(
    proc.time()["elapsed"]
  )
  
  set.seed(iteration_seed)
  
  sampled_cells <- unlist(
    lapply(
      names(cells_by_donor),
      function(donor_name) {
        
        donor_cells <- unlist(
          cells_by_donor[donor_name],
          use.names = FALSE
        )
        
        donor_target <- as.integer(
          target_cells_by_donor[donor_name]
        )
        
        sample(
          donor_cells,
          size = donor_target,
          replace = FALSE
        )
      }
    ),
    use.names = FALSE
  )
  
  if (length(sampled_cells) != matched_total_cells) {
    stop(
      "Iteration ", iteration_number,
      " sampled an incorrect number of cells."
    )
  }
  
  control_counts <- raw_counts[
    ,
    sampled_cells,
    drop = FALSE
  ]
  
  control_metadata <- metadata_source[
    sampled_cells,
    ,
    drop = FALSE
  ]
  
  control_metadata$donor_id_analysis <- as.character(
    control_metadata[, donor_column]
  )
  
  observed_donor_counts <- table(
    control_metadata$donor_id_analysis
  )
  
  expected_donor_counts <- target_cells_by_donor[
    names(observed_donor_counts)
  ]
  
  if (
    !identical(
      as.integer(observed_donor_counts),
      as.integer(expected_donor_counts)
    )
  ) {
    stop(
      "Iteration ", iteration_number,
      " did not reproduce the target donor counts."
    )
  }
  
  control_obj <- CreateSeuratObject(
    counts = control_counts,
    meta.data = control_metadata,
    assay = "RNA",
    project = paste0(
      "MatchedUnbalanced_",
      iteration_number
    ),
    min.cells = 0,
    min.features = 0
  )
  
  DefaultAssay(control_obj) <- "RNA"
  
  control_obj[["RNA"]] <- split(
    control_obj[["RNA"]],
    f = control_obj$donor_id_analysis
  )
  
  set.seed(iteration_seed)
  
  control_obj <- NormalizeData(
    object = control_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  
  control_obj <- FindVariableFeatures(
    object = control_obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  
  control_obj <- ScaleData(
    object = control_obj,
    features = VariableFeatures(control_obj),
    verbose = FALSE
  )
  
  control_obj <- RunPCA(
    object = control_obj,
    features = VariableFeatures(control_obj),
    npcs = 20,
    seed.use = iteration_seed,
    verbose = FALSE
  )
  
  set.seed(iteration_seed)
  
  control_obj <- IntegrateLayers(
    object = control_obj,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca",
    dims = 1:20,
    verbose = FALSE
  )
  
  control_obj <- FindNeighbors(
    object = control_obj,
    reduction = "integrated.rpca",
    dims = 1:20,
    k.param = 20,
    verbose = FALSE
  )
  
  set.seed(iteration_seed)
  
  control_obj <- FindClusters(
    object = control_obj,
    resolution = 0.6,
    algorithm = 3,
    random.seed = iteration_seed,
    verbose = FALSE
  )
  
  control_labels <- as.character(
    Idents(control_obj)
  )
  
  names(control_labels) <- colnames(control_obj)
  
  evaluation_cells <- sampled_cells[
    sampled_cells %in% names(control_labels) &
      sampled_cells %in% names(full_reference_vector)
  ]
  
  if (length(evaluation_cells) != matched_total_cells) {
    stop(
      "Iteration ", iteration_number,
      " could not align every sampled cell."
    )
  }
  
  reference_labels_subset <-
    full_reference_vector[evaluation_cells]
  
  control_labels_subset <-
    control_labels[evaluation_cells]
  
  ari_value <- mclust::adjustedRandIndex(
    reference_labels_subset,
    control_labels_subset
  )
  
  nmi_value <- aricode::NMI(
    reference_labels_subset,
    control_labels_subset
  )
  
  runtime_seconds <- unname(
    proc.time()["elapsed"]
  ) - iteration_start
  
  iteration_metrics <- data.frame(
    iteration = iteration_number,
    seed = iteration_seed,
    total_sampled_cells = matched_total_cells,
    control_clusters_recovered = length(
      unique(control_labels_subset)
    ),
    reference_clusters_present = length(
      unique(reference_labels_subset)
    ),
    ARI_vs_full_reference = round(ari_value, 4),
    NMI_vs_full_reference = round(nmi_value, 4),
    runtime_seconds = round(runtime_seconds, 2),
    status = "Completed",
    stringsAsFactors = FALSE
  )
  
  sampled_cell_record <- data.frame(
    iteration = iteration_number,
    seed = iteration_seed,
    cell_barcode = sampled_cells,
    donor_id = control_metadata$donor_id_analysis,
    reference_cluster =
      full_reference_vector[sampled_cells],
    stringsAsFactors = FALSE
  )
  
  rm(control_obj)
  invisible(gc())
  
  list(
    metrics = iteration_metrics,
    sampled_cells = sampled_cell_record
  )
}

cat("Matched-size unbalanced iteration function defined.\n")





# ============================================================
# 8. Run one pilot control iteration
# ============================================================

cat("Running pilot matched-size unbalanced iteration...\n")

pilot_control_result <-
  run_unbalanced_control_iteration(
    iteration_number = 1,
    iteration_seed = iteration_seeds[1]
  )

print(
  pilot_control_result$metrics,
  row.names = FALSE
)

cat("\nPilot sampled donor counts:\n")

print(
  table(
    pilot_control_result$sampled_cells$donor_id
  )
)

write.csv(
  pilot_control_result$metrics,
  file.path(
    output_dir,
    "pilot_matched_size_unbalanced_metrics.csv"
  ),
  row.names = FALSE
)

cat("Pilot control iteration completed successfully.\n")



# ============================================================
# 9. Run the remaining matched-size control iterations
# ============================================================

results_list <- vector(
  "list",
  number_of_iterations
)

sampled_cell_list <- vector(
  "list",
  number_of_iterations
)

failure_list <- list()

results_list[[1]] <- pilot_control_result$metrics
sampled_cell_list[[1]] <- pilot_control_result$sampled_cells

for (iteration_number in 2:number_of_iterations) {
  
  iteration_seed <- iteration_seeds[iteration_number]
  
  cat(
    "\nRunning control iteration",
    iteration_number,
    "of",
    number_of_iterations,
    "with seed",
    iteration_seed,
    "...\n"
  )
  
  iteration_output <- tryCatch(
    {
      run_unbalanced_control_iteration(
        iteration_number = iteration_number,
        iteration_seed = iteration_seed
      )
    },
    error = function(e) {
      
      failure_list[[length(failure_list) + 1]] <<-
        data.frame(
          iteration = iteration_number,
          seed = iteration_seed,
          error_message = conditionMessage(e),
          stringsAsFactors = FALSE
        )
      
      cat(
        "Iteration",
        iteration_number,
        "failed:",
        conditionMessage(e),
        "\n"
      )
      
      NULL
    }
  )
  
  if (!is.null(iteration_output)) {
    
    results_list[[iteration_number]] <-
      iteration_output$metrics
    
    sampled_cell_list[[iteration_number]] <-
      iteration_output$sampled_cells
    
    cat(
      "Completed iteration",
      iteration_number,
      ": ARI =",
      iteration_output$metrics$ARI_vs_full_reference,
      ", NMI =",
      iteration_output$metrics$NMI_vs_full_reference,
      ", clusters =",
      iteration_output$metrics$control_clusters_recovered,
      "\n"
    )
  }
  
  if (iteration_number %% 5 == 0) {
    
    checkpoint_results <- bind_rows(
      results_list[
        !vapply(
          results_list,
          is.null,
          logical(1)
        )
      ]
    )
    
    write.csv(
      checkpoint_results,
      file.path(
        output_dir,
        paste0(
          "checkpoint_metrics_after_iteration_",
          iteration_number,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    
    cat(
      "Checkpoint saved after iteration",
      iteration_number,
      "\n"
    )
  }
  
  invisible(gc())
}

completed_iteration_count <- sum(
  !vapply(
    results_list,
    is.null,
    logical(1)
  )
)

cat(
  "\nControl loop finished. Completed:",
  completed_iteration_count,
  "of",
  number_of_iterations,
  "iterations.\n"
)




# ============================================================
# 10. Combine and save completed control results
# ============================================================

completed_results <- bind_rows(
  results_list[
    !vapply(
      results_list,
      is.null,
      logical(1)
    )
  ]
)

completed_sampled_cells <- bind_rows(
  sampled_cell_list[
    !vapply(
      sampled_cell_list,
      is.null,
      logical(1)
    )
  ]
)

stopifnot(
  nrow(completed_results) == 50,
  nrow(completed_sampled_cells) ==
    50 * matched_total_cells,
  all(completed_results$status == "Completed")
)

write.csv(
  completed_results,
  file.path(
    output_dir,
    "Table_matched_size_unbalanced_all_iterations.csv"
  ),
  row.names = FALSE
)

write.csv(
  completed_sampled_cells,
  file.path(
    output_dir,
    "sampled_cells_all_unbalanced_iterations.csv"
  ),
  row.names = FALSE
)

if (length(failure_list) > 0) {
  
  failed_iterations <- bind_rows(
    failure_list
  )
  
  write.csv(
    failed_iterations,
    file.path(
      output_dir,
      "failed_iterations.csv"
    ),
    row.names = FALSE
  )
}

cat(
  "Saved", nrow(completed_results),
  "iteration records and",
  nrow(completed_sampled_cells),
  "sampled-cell records.\n"
)

cat(
  "Failed iterations:",
  length(failure_list),
  "\n"
)



# ============================================================
# 11. Locate the completed balanced-analysis results
# ============================================================

balanced_result_candidates <- c(
  file.path(
    "repeated_donor_balanced_sensitivity_outputs",
    "Table_repeated_balanced_sensitivity_all_iterations.csv"
  ),
  file.path(
    "repeated_donor_balanced_sensitivity_outputs(1)",
    "Table_repeated_balanced_sensitivity_all_iterations.csv"
  ),
  "Table_repeated_balanced_sensitivity_all_iterations.csv"
)

available_balanced_files <-
  balanced_result_candidates[
    file.exists(balanced_result_candidates)
  ]

if (length(available_balanced_files) == 0) {
  stop(
    "The balanced 50-iteration results file was not found."
  )
}

balanced_results_file <-
  available_balanced_files[1]

cat(
  "Balanced results file found:",
  balanced_results_file,
  "\n"
)





balanced_results <- read.csv(
  balanced_results_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_balanced_columns <- c(
  "iteration",
  "ARI_vs_full_reference",
  "NMI_vs_full_reference",
  "balanced_clusters_recovered"
)

missing_balanced_columns <- setdiff(
  required_balanced_columns,
  colnames(balanced_results)
)

if (length(missing_balanced_columns) > 0) {
  stop(
    "Missing balanced-result columns: ",
    paste(missing_balanced_columns, collapse = ", ")
  )
}

stopifnot(
  nrow(balanced_results) == 50
)

cat(
  "Balanced results verified:",
  nrow(balanced_results),
  "iterations.\n"
)

print(
  colnames(balanced_results)
)


# ============================================================
# 12. Create balanced-versus-unbalanced comparison
# ============================================================

balanced_standardized <- data.frame(
  iteration = balanced_results$iteration,
  sampling_design = "Donor-balanced",
  total_cells =
    balanced_results$total_balanced_cells,
  ARI =
    balanced_results$ARI_vs_full_reference,
  NMI =
    balanced_results$NMI_vs_full_reference,
  clusters_recovered =
    balanced_results$balanced_clusters_recovered
)

unbalanced_standardized <- data.frame(
  iteration = completed_results$iteration,
  sampling_design = "Matched-size unbalanced",
  total_cells =
    completed_results$total_sampled_cells,
  ARI =
    completed_results$ARI_vs_full_reference,
  NMI =
    completed_results$NMI_vs_full_reference,
  clusters_recovered =
    completed_results$control_clusters_recovered
)

comparison_all_iterations <- bind_rows(
  balanced_standardized,
  unbalanced_standardized
)

stopifnot(
  nrow(comparison_all_iterations) == 100,
  all(
    comparison_all_iterations$total_cells ==
      matched_total_cells
  )
)

comparison_summary <-
  comparison_all_iterations %>%
  group_by(sampling_design) %>%
  summarise(
    iterations = n(),
    mean_ARI = round(mean(ARI), 4),
    sd_ARI = round(sd(ARI), 4),
    median_ARI = round(median(ARI), 4),
    mean_NMI = round(mean(NMI), 4),
    sd_NMI = round(sd(NMI), 4),
    median_NMI = round(median(NMI), 4),
    mean_clusters = round(
      mean(clusters_recovered),
      2
    ),
    sd_clusters = round(
      sd(clusters_recovered),
      2
    ),
    median_clusters = median(
      clusters_recovered
    ),
    minimum_clusters = min(
      clusters_recovered
    ),
    maximum_clusters = max(
      clusters_recovered
    ),
    .groups = "drop"
  )

write.csv(
  comparison_all_iterations,
  file.path(
    output_dir,
    "Table_balanced_vs_unbalanced_all_iterations.csv"
  ),
  row.names = FALSE
)

write.csv(
  comparison_summary,
  file.path(
    output_dir,
    "Table_balanced_vs_unbalanced_summary.csv"
  ),
  row.names = FALSE
)

print(
  as.data.frame(comparison_summary),
  row.names = FALSE
)


# ============================================================
# 13. Correct statistical comparison
# ============================================================
cliffs_delta <- function(x, y) {
  
  pairwise_differences <- outer(
    x,
    y,
    FUN = "-"
  )
  
  (
    sum(pairwise_differences > 0) -
      sum(pairwise_differences < 0)
  ) / length(pairwise_differences)
}

compare_metric_corrected <- function(
    metric_name,
    balanced_values,
    unbalanced_values
) {
  
  welch_result <- t.test(
    balanced_values,
    unbalanced_values,
    var.equal = FALSE
  )
  
  confidence_interval <- unname(
    welch_result$conf.int
  )
  
  delta_value <- cliffs_delta(
    balanced_values,
    unbalanced_values
  )
  
  delta_magnitude <- cut(
    abs(delta_value),
    breaks = c(
      -Inf,
      0.147,
      0.33,
      0.474,
      Inf
    ),
    labels = c(
      "Negligible",
      "Small",
      "Medium",
      "Large"
    ),
    right = FALSE
  )
  
  data.frame(
    metric = metric_name,
    balanced_mean = mean(balanced_values),
    unbalanced_mean = mean(unbalanced_values),
    mean_difference =
      mean(balanced_values) -
      mean(unbalanced_values),
    Welch_CI95_lower = confidence_interval[1],
    Welch_CI95_upper = confidence_interval[2],
    Welch_p_value = welch_result$p.value,
    Wilcoxon_p_value = wilcox.test(
      balanced_values,
      unbalanced_values,
      exact = FALSE
    )$p.value,
    cliffs_delta = delta_value,
    effect_magnitude =
      as.character(delta_magnitude),
    stringsAsFactors = FALSE
  )
}

statistical_comparison_corrected <- bind_rows(
  compare_metric_corrected(
    "ARI",
    balanced_standardized$ARI,
    unbalanced_standardized$ARI
  ),
  compare_metric_corrected(
    "NMI",
    balanced_standardized$NMI,
    unbalanced_standardized$NMI
  ),
  compare_metric_corrected(
    "Clusters recovered",
    balanced_standardized$clusters_recovered,
    unbalanced_standardized$clusters_recovered
  )
)

statistical_comparison_corrected$
  BH_adjusted_Wilcoxon_p_value <- p.adjust(
    statistical_comparison_corrected$
      Wilcoxon_p_value,
    method = "BH"
  )

statistical_comparison_corrected$
  BH_adjusted_Welch_p_value <- p.adjust(
    statistical_comparison_corrected$
      Welch_p_value,
    method = "BH"
  )

write.csv(
  statistical_comparison_corrected,
  file.path(
    output_dir,
    "Table_balanced_vs_unbalanced_statistical_comparison_corrected.csv"
  ),
  row.names = FALSE
)

print(
  statistical_comparison_corrected,
  row.names = FALSE,
  digits = 6
)



# ============================================================
# 14. Generate balanced-versus-unbalanced comparison figure
# ============================================================

comparison_long <- comparison_all_iterations %>%
  pivot_longer(
    cols = c(
      ARI,
      NMI,
      clusters_recovered
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    sampling_label = factor(
      sampling_design,
      levels = c(
        "Donor-balanced",
        "Matched-size unbalanced"
      ),
      labels = c(
        "Balanced",
        "Unbalanced"
      )
    ),
    metric = factor(
      metric,
      levels = c(
        "ARI",
        "NMI",
        "clusters_recovered"
      ),
      labels = c(
        "ARI",
        "NMI",
        "Clusters recovered"
      )
    )
  )

comparison_figure <- ggplot(
  comparison_long,
  aes(
    x = sampling_label,
    y = value,
    fill = sampling_design
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA,
    alpha = 0.85
  ) +
  geom_jitter(
    width = 0.12,
    size = 1.5,
    alpha = 0.65
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    nrow = 1
  ) +
  scale_fill_manual(
    values = c(
      "Donor-balanced" = "#0072B2",
      "Matched-size unbalanced" = "#D55E00"
    )
  ) +
  labs(
    title = "Balanced versus unbalanced donor sampling",
    subtitle = paste0(
      "Fifty RPCA reruns with ",
      matched_total_cells,
      " cells per iteration"
    ),
    x = NULL,
    y = "Metric value"
  ) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    strip.text = element_text(
      face = "bold"
    ),
    axis.text.x = element_text(
      angle = 0,
      hjust = 0.5
    )
  )

print(comparison_figure)

ggsave(
  file.path(
    output_dir,
    "Figure_balanced_vs_matched_unbalanced_comparison.png"
  ),
  comparison_figure,
  width = 12,
  height = 4.8,
  dpi = 400,
  bg = "white"
)

ggsave(
  file.path(
    output_dir,
    "Figure_balanced_vs_matched_unbalanced_comparison.pdf"
  ),
  comparison_figure,
  width = 12,
  height = 4.8,
  bg = "white"
)

cat("Final comparison figure saved.\n")




# ============================================================
# 15. Save analysis summary and session information
# ============================================================

summary_file <- file.path(
  output_dir,
  "matched_size_unbalanced_control_summary.txt"
)

capture.output(
  {
    cat(
      "MATCHED-SIZE UNBALANCED DONOR CONTROL\n",
      "=======================================\n\n"
    )
    
    cat(
      "Iterations per sampling design: 50\n"
    )
    
    cat(
      "Cells per iteration:",
      matched_total_cells,
      "\n\n"
    )
    
    cat("UNBALANCED SAMPLING DESIGN\n")
    print(sampling_design, row.names = FALSE)
    
    cat("\nDESCRIPTIVE COMPARISON\n")
    print(
      as.data.frame(comparison_summary),
      row.names = FALSE
    )
    
    cat("\nSTATISTICAL COMPARISON\n")
    print(
      statistical_comparison_corrected,
      row.names = FALSE
    )
    
    cat(
      "\nINTERPRETATION\n",
      "Both sampling designs contained exactly ",
      matched_total_cells,
      " cells per iteration. Donor-balanced sampling ",
      "produced higher ARI and NMI agreement with the ",
      "full-data reference and recovered more clusters ",
      "than proportional unbalanced sampling. Therefore, ",
      "the balanced-analysis results cannot be explained ",
      "only by reducing the total number of cells; donor ",
      "representation also affected clustering preservation.\n",
      sep = ""
    )
  },
  file = summary_file
)

session_info_file <- file.path(
  output_dir,
  "matched_size_unbalanced_session_info.txt"
)

capture.output(
  sessionInfo(),
  file = session_info_file
)

cat(
  "Summary saved:",
  file.exists(summary_file),
  "\nSession information saved:",
  file.exists(session_info_file),
  "\n"
)





invalid_statistical_file <- file.path(
  output_dir,
  "Table_balanced_vs_unbalanced_statistical_comparison.csv"
)

corrected_statistical_file <- file.path(
  output_dir,
  "Table_balanced_vs_unbalanced_statistical_comparison_corrected.csv"
)

stopifnot(
  file.exists(corrected_statistical_file)
)

if (file.exists(invalid_statistical_file)) {
  file.remove(invalid_statistical_file)
}

cat(
  "Invalid preliminary table removed:",
  !file.exists(invalid_statistical_file),
  "\nCorrected table retained:",
  file.exists(corrected_statistical_file),
  "\n"
)