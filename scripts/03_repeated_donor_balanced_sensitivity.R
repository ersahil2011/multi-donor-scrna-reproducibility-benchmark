# ============================================================
# 03_repeated_donor_balanced_sensitivity.R
#
# Purpose:
# Perform 50 repeated donor-balanced sensitivity analyses.
# For each iteration, equal numbers of cells are sampled from
# each donor, the RPCA-based clustering workflow is rebuilt,
# and agreement with the retained full-data clustering solution
# at resolution 0.6 is quantified using ARI and NMI.
#
# Original RDS objects are read but never modified.
# ============================================================


# ============================================================
# 1. Load required packages
# ============================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "mclust",
  "aricode",
  "dplyr",
  "ggplot2",
  "future"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(Matrix)
  library(mclust)
  library(aricode)
  library(dplyr)
  library(ggplot2)
  library(future)
})

# Use sequential processing to avoid high memory usage.
future::plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

cat("Required packages loaded successfully.\n")


# ============================================================
# 2. Create clean output folder
# ============================================================

output_dir <- "repeated_donor_balanced_sensitivity_outputs"

if (dir.exists(output_dir)) {
  unlink(output_dir, recursive = TRUE)
}

dir.create(output_dir)

cat("Fresh output folder created:", output_dir, "\n")


# ============================================================
# 3. Locate and load required RDS objects
# ============================================================

post_filter_candidates <- c(
  "sce_after_doublet_removal.rds",
  "sce_after_doublet_removal(1).rds"
)

reference_candidates <- c(
  "sce_final_clusters_locked_rpca_pcs20_res06.rds",
  "sce_final_clusters_locked_rpca_pcs20_res06(1).rds"
)

post_filter_file <- post_filter_candidates[
  file.exists(post_filter_candidates)
][1]

reference_file <- reference_candidates[
  file.exists(reference_candidates)
][1]

if (is.na(post_filter_file)) {
  stop("Post-doublet-removal object not found in the project folder.")
}

if (is.na(reference_file)) {
  stop("Final reference clustering object not found in the project folder.")
}

cat("Post-filtering object used:", post_filter_file, "\n")
cat("Reference clustering object used:", reference_file, "\n")

post_obj <- readRDS(post_filter_file)
reference_obj <- readRDS(reference_file)

if (!inherits(post_obj, "Seurat")) {
  stop("The post-doublet-removal object is not a Seurat object.")
}

if (!inherits(reference_obj, "Seurat")) {
  stop("The final reference clustering object is not a Seurat object.")
}

DefaultAssay(post_obj) <- "RNA"

cat("\nPost-filtering object dimensions:\n")
print(dim(post_obj))

cat("\nReference object dimensions:\n")
print(dim(reference_obj))


# ============================================================
# 4. Detect donor metadata column
# ============================================================

detect_donor_column <- function(metadata) {
  
  preferred_columns <- c(
    "sample_id",
    "donor_id",
    "donor",
    "orig.ident",
    "sample"
  )
  
  detected <- preferred_columns[
    preferred_columns %in% colnames(metadata)
  ]
  
  if (length(detected) > 0) {
    return(detected[1])
  }
  
  fallback <- grep(
    pattern = "sample|donor|orig.ident",
    x = colnames(metadata),
    value = TRUE,
    ignore.case = TRUE
  )
  
  if (length(fallback) > 0) {
    return(fallback[1])
  }
  
  return(NA_character_)
}

donor_column <- detect_donor_column(post_obj@meta.data)

if (is.na(donor_column)) {
  stop("No donor metadata column could be identified.")
}

cat("\nDetected donor metadata column:", donor_column, "\n")

original_donor_counts <- table(post_obj@meta.data[[donor_column]])

cat("\nFinal retained cell counts by donor:\n")
print(original_donor_counts)

if (length(original_donor_counts) != 5) {
  stop("The detected donor column does not contain exactly five donors.")
}

write.csv(
  data.frame(
    donor = names(original_donor_counts),
    retained_cells = as.integer(original_donor_counts)
  ),
  file.path(output_dir, "donor_counts_full_final_object.csv"),
  row.names = FALSE
)


# ============================================================
# 5. Recover and align donor-specific raw count layers
# ============================================================

rna_layer_names <- tryCatch(
  {
    Layers(post_obj[["RNA"]])
  },
  error = function(e) {
    character(0)
  }
)

count_layer_names <- rna_layer_names[
  grepl("^counts($|\\.)", rna_layer_names, ignore.case = TRUE)
]

cat("\nCounts-related RNA layers detected:\n")
print(count_layer_names)

if (length(count_layer_names) == 0) {
  stop("No raw count layers were detected in the RNA assay.")
}

feature_universe <- rownames(post_obj)

cat("\nCommon analytical feature universe:", length(feature_universe), "features\n")

extract_and_align_count_layer <- function(layer_name,
                                          feature_universe,
                                          seurat_object) {
  
  cat("Extracting and aligning layer:", layer_name, "\n")
  
  layer_matrix <- LayerData(
    object = seurat_object,
    assay = "RNA",
    layer = layer_name
  )
  
  unexpected_features <- setdiff(
    rownames(layer_matrix),
    feature_universe
  )
  
  if (length(unexpected_features) > 0) {
    stop(
      paste0(
        "Layer ", layer_name,
        " contains features outside the common feature universe."
      )
    )
  }
  
  missing_features <- setdiff(
    feature_universe,
    rownames(layer_matrix)
  )
  
  if (length(missing_features) > 0) {
    
    zero_matrix <- Matrix::Matrix(
      data = 0,
      nrow = length(missing_features),
      ncol = ncol(layer_matrix),
      sparse = TRUE,
      dimnames = list(
        missing_features,
        colnames(layer_matrix)
      )
    )
    
    layer_matrix <- rbind(layer_matrix, zero_matrix)
  }
  
  layer_matrix <- layer_matrix[
    feature_universe,
    ,
    drop = FALSE
  ]
  
  cat(
    "Aligned layer dimensions:",
    paste(dim(layer_matrix), collapse = " x "),
    "\n"
  )
  
  return(layer_matrix)
}

count_matrix_list <- lapply(
  count_layer_names,
  extract_and_align_count_layer,
  feature_universe = feature_universe,
  seurat_object = post_obj
)

raw_counts <- do.call(cbind, count_matrix_list)

cat("\nCombined aligned raw count matrix dimensions:\n")
print(dim(raw_counts))

if (anyDuplicated(colnames(raw_counts)) > 0) {
  stop("Duplicate cell barcodes were detected after combining count layers.")
}

common_cells <- intersect(
  colnames(raw_counts),
  rownames(post_obj@meta.data)
)

if (length(common_cells) == 0) {
  stop("No common cells were found between raw counts and metadata.")
}

raw_counts <- raw_counts[, common_cells, drop = FALSE]
metadata_source <- post_obj@meta.data[common_cells, , drop = FALSE]

aligned_donor_counts <- table(metadata_source[[donor_column]])

cat("\nDonor distribution after count-metadata alignment:\n")
print(aligned_donor_counts)

if (length(aligned_donor_counts) != 5) {
  stop("The aligned dataset does not contain exactly five donors.")
}

cells_per_donor <- min(as.integer(aligned_donor_counts))
total_balanced_cells <- cells_per_donor * length(aligned_donor_counts)

cat("\nRepeated balanced-sampling design:\n")
cat("Cells per donor per iteration:", cells_per_donor, "\n")
cat("Total balanced cells per iteration:", total_balanced_cells, "\n")


# ============================================================
# 6. Identify full-data reference cluster assignments
# ============================================================

reference_metadata <- reference_obj@meta.data

reference_column_candidates <- c(
  "cluster_final",
  "RNA_snn_res.0.6",
  "integrated_snn_res.0.6"
)

reference_cluster_column <- reference_column_candidates[
  reference_column_candidates %in% colnames(reference_metadata)
][1]

if (is.na(reference_cluster_column)) {
  
  fallback_columns <- grep(
    pattern = "cluster_final|res\\.0\\.6$|res\\.0.6$",
    x = colnames(reference_metadata),
    value = TRUE,
    ignore.case = TRUE
  )
  
  if (length(fallback_columns) == 0) {
    stop("The full-data resolution-0.6 reference clustering column was not detected.")
  }
  
  reference_cluster_column <- fallback_columns[1]
}

cat("\nFull-data reference cluster column:", reference_cluster_column, "\n")

full_reference_vector <- as.character(
  reference_metadata[[reference_cluster_column]]
)

names(full_reference_vector) <- rownames(reference_metadata)


# ============================================================
# 7. Set repeated sensitivity-analysis parameters
# ============================================================

number_of_iterations <- 50
iteration_seeds <- 1001:(1000 + number_of_iterations)

write.csv(
  data.frame(
    iteration = seq_len(number_of_iterations),
    seed = iteration_seeds
  ),
  file.path(output_dir, "iteration_seed_record.csv"),
  row.names = FALSE
)

cells_by_donor <- split(
  rownames(metadata_source),
  metadata_source[[donor_column]]
)

results_list <- vector("list", number_of_iterations)
sampled_cell_list <- vector("list", number_of_iterations)
failure_list <- list()


# ============================================================
# 8. Define one repeated balanced RPCA analysis function
# ============================================================

run_balanced_iteration <- function(iteration_number, iteration_seed) {
  
  set.seed(iteration_seed)
  
  sampled_cells <- unlist(
    lapply(
      cells_by_donor,
      function(cell_vector) {
        sample(
          x = cell_vector,
          size = cells_per_donor,
          replace = FALSE
        )
      }
    ),
    use.names = FALSE
  )
  
  balanced_counts <- raw_counts[, sampled_cells, drop = FALSE]
  balanced_metadata <- metadata_source[sampled_cells, , drop = FALSE]
  balanced_metadata$donor_id_analysis <- as.character(
    balanced_metadata[[donor_column]]
  )
  
  balanced_obj <- CreateSeuratObject(
    counts = balanced_counts,
    meta.data = balanced_metadata,
    assay = "RNA",
    project = paste0("BalancedIteration_", iteration_number),
    min.cells = 0,
    min.features = 0
  )
  
  DefaultAssay(balanced_obj) <- "RNA"
  
  balanced_obj[["RNA"]] <- split(
    balanced_obj[["RNA"]],
    f = balanced_obj$donor_id_analysis
  )
  
  set.seed(iteration_seed)
  
  balanced_obj <- NormalizeData(
    object = balanced_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  
  balanced_obj <- FindVariableFeatures(
    object = balanced_obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  
  balanced_obj <- ScaleData(
    object = balanced_obj,
    features = VariableFeatures(balanced_obj),
    verbose = FALSE
  )
  
  balanced_obj <- RunPCA(
    object = balanced_obj,
    features = VariableFeatures(balanced_obj),
    npcs = 20,
    seed.use = iteration_seed,
    verbose = FALSE
  )
  
  set.seed(iteration_seed)
  
  balanced_obj <- IntegrateLayers(
    object = balanced_obj,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca",
    dims = 1:20,
    verbose = FALSE
  )
  
  balanced_obj <- FindNeighbors(
    object = balanced_obj,
    reduction = "integrated.rpca",
    dims = 1:20,
    k.param = 20,
    verbose = FALSE
  )
  
  set.seed(iteration_seed)
  
  balanced_obj <- FindClusters(
    object = balanced_obj,
    resolution = 0.6,
    algorithm = 3,
    random.seed = iteration_seed,
    verbose = FALSE
  )
  
  balanced_labels <- as.character(Idents(balanced_obj))
  names(balanced_labels) <- colnames(balanced_obj)
  
  evaluation_cells <- intersect(
    names(balanced_labels),
    names(full_reference_vector)
  )
  
  if (length(evaluation_cells) != total_balanced_cells) {
    stop(
      paste0(
        "Iteration ", iteration_number,
        " did not retain all balanced cells for reference comparison."
      )
    )
  }
  
  reference_labels_subset <- full_reference_vector[evaluation_cells]
  balanced_labels_subset <- balanced_labels[evaluation_cells]
  
  ari_value <- mclust::adjustedRandIndex(
    reference_labels_subset,
    balanced_labels_subset
  )
  
  nmi_value <- aricode::NMI(
    reference_labels_subset,
    balanced_labels_subset
  )
  
  number_of_balanced_clusters <- length(
    unique(balanced_labels_subset)
  )
  
  number_of_reference_clusters_present <- length(
    unique(reference_labels_subset)
  )
  
  iteration_result <- data.frame(
    iteration = iteration_number,
    seed = iteration_seed,
    cells_per_donor = cells_per_donor,
    total_balanced_cells = total_balanced_cells,
    balanced_clusters_recovered = number_of_balanced_clusters,
    reference_clusters_present_in_balanced_cells =
      number_of_reference_clusters_present,
    ARI_vs_full_reference = round(ari_value, 4),
    NMI_vs_full_reference = round(nmi_value, 4),
    status = "Completed",
    stringsAsFactors = FALSE
  )
  
  sampled_cell_record <- data.frame(
    iteration = iteration_number,
    seed = iteration_seed,
    cell_barcode = sampled_cells,
    donor_id = balanced_metadata$donor_id_analysis,
    stringsAsFactors = FALSE
  )
  
  rm(balanced_obj)
  invisible(gc())
  
  return(
    list(
      metrics = iteration_result,
      sampled_cells = sampled_cell_record
    )
  )
}


# ============================================================
# 9. Execute 50 repeated balanced sensitivity iterations
# ============================================================

analysis_start_time <- Sys.time()

for (iteration_number in seq_len(number_of_iterations)) {
  
  iteration_seed <- iteration_seeds[iteration_number]
  
  cat(
    "\nRunning balanced sensitivity iteration",
    iteration_number,
    "of",
    number_of_iterations,
    "with seed",
    iteration_seed,
    "...\n"
  )
  
  iteration_output <- tryCatch(
    {
      run_balanced_iteration(
        iteration_number = iteration_number,
        iteration_seed = iteration_seed
      )
    },
    error = function(e) {
      
      failure_list[[length(failure_list) + 1]] <<- data.frame(
        iteration = iteration_number,
        seed = iteration_seed,
        error_message = conditionMessage(e),
        stringsAsFactors = FALSE
      )
      
      return(NULL)
    }
  )
  
  if (!is.null(iteration_output)) {
    
    results_list[[iteration_number]] <- iteration_output$metrics
    sampled_cell_list[[iteration_number]] <- iteration_output$sampled_cells
    
    cat(
      "Completed iteration",
      iteration_number,
      ": ARI =",
      iteration_output$metrics$ARI_vs_full_reference,
      ", NMI =",
      iteration_output$metrics$NMI_vs_full_reference,
      ", clusters =",
      iteration_output$metrics$balanced_clusters_recovered,
      "\n"
    )
  }
  
  # Save checkpoint after every five iterations.
  if (iteration_number %% 5 == 0) {
    
    completed_results_checkpoint <- bind_rows(
      results_list[!vapply(results_list, is.null, logical(1))]
    )
    
    write.csv(
      completed_results_checkpoint,
      file.path(
        output_dir,
        paste0("checkpoint_metrics_after_iteration_", iteration_number, ".csv")
      ),
      row.names = FALSE
    )
    
    cat("Checkpoint saved after iteration", iteration_number, "\n")
  }
}

analysis_end_time <- Sys.time()
total_runtime_minutes <- as.numeric(
  difftime(
    analysis_end_time,
    analysis_start_time,
    units = "mins"
  )
)


# ============================================================
# 10. Combine and export final iteration-level results
# ============================================================

completed_results <- bind_rows(
  results_list[!vapply(results_list, is.null, logical(1))]
)

completed_sampled_cells <- bind_rows(
  sampled_cell_list[!vapply(sampled_cell_list, is.null, logical(1))]
)

if (nrow(completed_results) == 0) {
  stop("No repeated sensitivity-analysis iterations completed successfully.")
}

write.csv(
  completed_results,
  file.path(output_dir, "Table_repeated_balanced_sensitivity_all_iterations.csv"),
  row.names = FALSE
)

write.csv(
  completed_sampled_cells,
  file.path(output_dir, "sampled_cells_all_iterations.csv"),
  row.names = FALSE
)

if (length(failure_list) > 0) {
  
  failed_iterations <- bind_rows(failure_list)
  
  write.csv(
    failed_iterations,
    file.path(output_dir, "failed_iterations.csv"),
    row.names = FALSE
  )
  
} else {
  
  failed_iterations <- data.frame()
}


# ============================================================
# 11. Create statistical summary table
# ============================================================

metric_summary <- data.frame(
  Metric = c(
    "ARI versus full-data reference",
    "NMI versus full-data reference",
    "Balanced clusters recovered"
  ),
  Mean = c(
    mean(completed_results$ARI_vs_full_reference),
    mean(completed_results$NMI_vs_full_reference),
    mean(completed_results$balanced_clusters_recovered)
  ),
  SD = c(
    sd(completed_results$ARI_vs_full_reference),
    sd(completed_results$NMI_vs_full_reference),
    sd(completed_results$balanced_clusters_recovered)
  ),
  Median = c(
    median(completed_results$ARI_vs_full_reference),
    median(completed_results$NMI_vs_full_reference),
    median(completed_results$balanced_clusters_recovered)
  ),
  Q1 = c(
    quantile(completed_results$ARI_vs_full_reference, 0.25),
    quantile(completed_results$NMI_vs_full_reference, 0.25),
    quantile(completed_results$balanced_clusters_recovered, 0.25)
  ),
  Q3 = c(
    quantile(completed_results$ARI_vs_full_reference, 0.75),
    quantile(completed_results$NMI_vs_full_reference, 0.75),
    quantile(completed_results$balanced_clusters_recovered, 0.75)
  ),
  Minimum = c(
    min(completed_results$ARI_vs_full_reference),
    min(completed_results$NMI_vs_full_reference),
    min(completed_results$balanced_clusters_recovered)
  ),
  Maximum = c(
    max(completed_results$ARI_vs_full_reference),
    max(completed_results$NMI_vs_full_reference),
    max(completed_results$balanced_clusters_recovered)
  )
)

numeric_columns <- c("Mean", "SD", "Median", "Q1", "Q3", "Minimum", "Maximum")
metric_summary[numeric_columns] <- lapply(
  metric_summary[numeric_columns],
  function(x) round(x, 4)
)

write.csv(
  metric_summary,
  file.path(output_dir, "Table_repeated_balanced_sensitivity_summary.csv"),
  row.names = FALSE
)

cluster_frequency <- completed_results %>%
  count(balanced_clusters_recovered, name = "number_of_iterations") %>%
  mutate(
    percentage_of_completed_iterations =
      round(100 * number_of_iterations / nrow(completed_results), 2)
  )

write.csv(
  cluster_frequency,
  file.path(output_dir, "Table_cluster_count_frequency_balanced_runs.csv"),
  row.names = FALSE
)

same_cluster_count_proportion <- mean(
  completed_results$balanced_clusters_recovered == 18
)


# ============================================================
# 12. Generate quantitative figures
# ============================================================

metric_long <- bind_rows(
  data.frame(
    iteration = completed_results$iteration,
    metric = "ARI",
    value = completed_results$ARI_vs_full_reference
  ),
  data.frame(
    iteration = completed_results$iteration,
    metric = "NMI",
    value = completed_results$NMI_vs_full_reference
  )
)

agreement_plot <- ggplot(
  metric_long,
  aes(x = iteration, y = value, linetype = metric, shape = metric)
) +
  geom_line(linewidth = 0.6) +
  geom_point(size = 2) +
  scale_y_continuous(limits = c(0, 1)) +
  labs(
    x = "Balanced sensitivity iteration",
    y = "Agreement with full-data reference clustering",
    linetype = "Metric",
    shape = "Metric"
  ) +
  theme_classic(base_size = 12)

ggsave(
  filename = file.path(output_dir, "Figure_ARI_NMI_across_balanced_iterations.png"),
  plot = agreement_plot,
  width = 8,
  height = 4.8,
  dpi = 600
)

cluster_frequency_plot <- ggplot(
  cluster_frequency,
  aes(x = factor(balanced_clusters_recovered),
      y = number_of_iterations)
) +
  geom_col() +
  labs(
    x = "Number of clusters recovered in balanced rerun",
    y = "Number of iterations"
  ) +
  theme_classic(base_size = 12)

ggsave(
  filename = file.path(output_dir, "Figure_cluster_count_frequency_balanced_runs.png"),
  plot = cluster_frequency_plot,
  width = 7,
  height = 4.8,
  dpi = 600
)


# ============================================================
# 13. Create a complete textual summary
# ============================================================

summary_text <- capture.output({
  
  cat("REPEATED DONOR-BALANCED SENSITIVITY ANALYSIS\n")
  cat("============================================\n\n")
  
  cat("Input post-filtering object:", post_filter_file, "\n")
  cat("Full-data reference object:", reference_file, "\n")
  cat("Donor metadata column:", donor_column, "\n")
  cat("Number of planned iterations:", number_of_iterations, "\n")
  cat("Number of completed iterations:", nrow(completed_results), "\n")
  cat("Number of failed iterations:", length(failure_list), "\n")
  cat("Cells sampled per donor per iteration:", cells_per_donor, "\n")
  cat("Total balanced cells per iteration:", total_balanced_cells, "\n")
  cat("Reference resolution: 0.6\n")
  cat("Total runtime in minutes:", round(total_runtime_minutes, 2), "\n\n")
  
  cat("FULL-DATA DONOR DISTRIBUTION\n")
  cat("----------------------------\n")
  print(aligned_donor_counts)
  
  cat("\nMETRIC SUMMARY ACROSS COMPLETED ITERATIONS\n")
  cat("------------------------------------------\n")
  print(metric_summary)
  
  cat("\nBALANCED CLUSTER-COUNT FREQUENCY\n")
  cat("--------------------------------\n")
  print(cluster_frequency)
  
  cat("\nProportion of completed balanced iterations recovering exactly 18 clusters:",
      round(same_cluster_count_proportion, 4), "\n")
  
  if (length(failure_list) > 0) {
    cat("\nFAILED ITERATIONS\n")
    cat("-----------------\n")
    print(failed_iterations)
  }
  
  cat("\nInterpretation note:\n")
  cat(
    "Each iteration independently balanced donor-wise cell representation ",
    "by sampling the same number of cells from each donor and rebuilding ",
    "the RPCA-based clustering workflow. ARI and NMI quantify agreement ",
    "between the balanced-rerun cluster assignment and the retained ",
    "full-data resolution-0.6 cluster assignment for the same cells. ",
    "Results should be interpreted as a sensitivity assessment of ",
    "clustering structure to donor representation, not as an independent ",
    "biological validation of the annotations.\n",
    sep = ""
  )
})

writeLines(
  summary_text,
  file.path(output_dir, "repeated_balanced_sensitivity_summary.txt")
)

session_text <- capture.output({
  cat("R SESSION INFORMATION\n")
  cat("=====================\n\n")
  print(sessionInfo())
})

writeLines(
  session_text,
  file.path(output_dir, "repeated_sensitivity_session_information.txt")
)


# ============================================================
# 14. Completion message
# ============================================================

cat("\n====================================================\n")
cat("REPEATED DONOR-BALANCED SENSITIVITY ANALYSIS COMPLETED\n")
cat("====================================================\n")
cat("Completed iterations:", nrow(completed_results), "of", number_of_iterations, "\n")
cat("Total runtime in minutes:", round(total_runtime_minutes, 2), "\n")
cat("Output folder created:", output_dir, "\n")
cat("Compress this folder as repeated_donor_balanced_sensitivity_outputs.zip\n")
cat("and upload it here for interpretation and manuscript integration.\n")