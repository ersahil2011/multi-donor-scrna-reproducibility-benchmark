# ============================================================
# 12_parameter_sensitivity.R
#
# Purpose:
# Evaluate clustering stability when varying:
# 1. Number of PCs: 10, 20, 30
# 2. Number of HVGs: 1000, 2000, 3000
# 3. Neighbour k: 10, 20, 30
#
# Parameters are varied one at a time around the retained
# setting: 2000 HVGs, 20 PCs, k = 20, resolution = 0.6.
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

cat("Parameter-sensitivity packages loaded successfully.\n")



# ============================================================
# 2. Load source objects
# ============================================================

output_dir <- file.path(
  "reviewer_revision_outputs",
  "12_parameter_sensitivity"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

input_files <- c(
  post_filter_object =
    "sce_after_doublet_removal.rds",
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

cat("Loading post-filtering object...\n")

post_obj <- readRDS(
  unname(
    input_files["post_filter_object"]
  )
)

cat("Loading reference object...\n")

reference_obj <- readRDS(
  unname(
    input_files["reference_object"]
  )
)

stopifnot(
  inherits(post_obj, "Seurat"),
  inherits(reference_obj, "Seurat"),
  identical(dim(post_obj), dim(reference_obj)),
  identical(colnames(post_obj), colnames(reference_obj)),
  "integrated.rpca" %in% Reductions(reference_obj)
)

DefaultAssay(post_obj) <- "RNA"
DefaultAssay(reference_obj) <- "RNA"

cat(
  "Objects loaded and aligned:",
  nrow(reference_obj), "genes and",
  ncol(reference_obj), "cells.\n"
)

cat("Output folder:", output_dir, "\n")




# ============================================================
# 3. Define reference clustering and sensitivity design
# ============================================================

reference_cluster_data <- FetchData(
  object = reference_obj,
  vars = "cluster_final",
  clean = "none"
)

reference_labels <- as.character(
  reference_cluster_data[, 1]
)

names(reference_labels) <- rownames(
  reference_cluster_data
)

reference_cluster_count <- length(
  unique(reference_labels)
)

available_rpca_dimensions <- ncol(
  Embeddings(
    reference_obj,
    reduction = "integrated.rpca"
  )
)

stopifnot(
  !anyNA(reference_labels),
  reference_cluster_count == 18,
  available_rpca_dimensions >= 30
)

baseline_HVGs <- 2000
baseline_PCs <- 20
baseline_k <- 20
fixed_resolution <- 0.6
clustering_seed <- 0

sensitivity_design <- data.frame(
  parameter = c(
    rep("PCs", 3),
    rep("HVGs", 3),
    rep("Neighbour_k", 3)
  ),
  value = c(
    10, 20, 30,
    1000, 2000, 3000,
    10, 20, 30
  ),
  baseline_value = c(
    rep(baseline_PCs, 3),
    rep(baseline_HVGs, 3),
    rep(baseline_k, 3)
  ),
  is_retained_setting = c(
    FALSE, TRUE, FALSE,
    FALSE, TRUE, FALSE,
    FALSE, TRUE, FALSE
  )
)

write.csv(
  sensitivity_design,
  file.path(
    output_dir,
    "Table_parameter_sensitivity_design.csv"
  ),
  row.names = FALSE
)

cat(
  "Reference clustering verified:",
  reference_cluster_count,
  "clusters.\n"
)

cat(
  "Available integrated RPCA dimensions:",
  available_rpca_dimensions,
  "\n"
)

print




# ============================================================
# 4. Define PC and neighbour-k sensitivity function
# ============================================================

run_graph_parameter_test <- function(
    parameter_name,
    parameter_value,
    dimensions_to_use,
    neighbour_k
) {
  
  test_start <- unname(
    proc.time()["elapsed"]
  )
  
  test_obj <- reference_obj
  
  graph_prefix <- paste0(
    "sensitivity_",
    tolower(parameter_name),
    "_",
    parameter_value
  )
  
  graph_names <- c(
    paste0(graph_prefix, "_nn"),
    paste0(graph_prefix, "_snn")
  )
  
  set.seed(clustering_seed)
  
  test_obj <- FindNeighbors(
    object = test_obj,
    reduction = "integrated.rpca",
    dims = dimensions_to_use,
    k.param = neighbour_k,
    graph.name = graph_names,
    verbose = FALSE
  )
  
  set.seed(clustering_seed)
  
  test_obj <- FindClusters(
    object = test_obj,
    graph.name = graph_names[2],
    resolution = fixed_resolution,
    algorithm = 3,
    random.seed = clustering_seed,
    verbose = FALSE
  )
  
  test_labels <- as.character(
    Idents(test_obj)
  )
  
  names(test_labels) <- colnames(test_obj)
  
  aligned_test_labels <- test_labels[
    names(reference_labels)
  ]
  
  if (anyNA(aligned_test_labels)) {
    stop(
      "Cell alignment failed for ",
      parameter_name,
      " = ",
      parameter_value
    )
  }
  
  runtime_seconds <- unname(
    proc.time()["elapsed"]
  ) - test_start
  
  result <- data.frame(
    parameter = parameter_name,
    value = parameter_value,
    HVGs = baseline_HVGs,
    PCs = length(dimensions_to_use),
    neighbour_k = neighbour_k,
    clusters_recovered = length(
      unique(aligned_test_labels)
    ),
    ARI_vs_retained_solution =
      mclust::adjustedRandIndex(
        reference_labels,
        aligned_test_labels
      ),
    NMI_vs_retained_solution =
      aricode::NMI(
        reference_labels,
        aligned_test_labels
      ),
    runtime_seconds = runtime_seconds,
    stringsAsFactors = FALSE
  )
  
  rm(test_obj)
  invisible(gc())
  
  result
}

cat("PC and neighbour-k sensitivity function defined.\n")



# ============================================================
# 5. Run PC and neighbour-k sensitivity tests
# ============================================================

pc_values <- c(10, 20, 30)
k_values <- c(10, 20, 30)

pc_results <- bind_rows(
  lapply(
    pc_values,
    function(pc_value) {
      
      cat(
        "Testing PC count:",
        pc_value,
        "\n"
      )
      
      run_graph_parameter_test(
        parameter_name = "PCs",
        parameter_value = pc_value,
        dimensions_to_use =
          seq_len(pc_value),
        neighbour_k = baseline_k
      )
    }
  )
)

k_results <- bind_rows(
  lapply(
    k_values,
    function(k_value) {
      
      cat(
        "Testing neighbour k:",
        k_value,
        "\n"
      )
      
      run_graph_parameter_test(
        parameter_name = "Neighbour_k",
        parameter_value = k_value,
        dimensions_to_use =
          seq_len(baseline_PCs),
        neighbour_k = k_value
      )
    }
  )
)

graph_sensitivity_results <- bind_rows(
  pc_results,
  k_results
)

write.csv(
  graph_sensitivity_results,
  file.path(
    output_dir,
    "Table_PC_and_neighbour_k_sensitivity.csv"
  ),
  row.names = FALSE
)

print(
  graph_sensitivity_results %>%
    mutate(
      ARI_vs_retained_solution = round(
        ARI_vs_retained_solution,
        4
      ),
      NMI_vs_retained_solution = round(
        NMI_vs_retained_solution,
        4
      ),
      runtime_seconds = round(
        runtime_seconds,
        2
      )
    ),
  row.names = FALSE
)

# ============================================================
# 6. Define raw-count layer alignment
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

if (length(count_layer_names) != 5) {
  stop(
    "Expected five count layers; detected ",
    length(count_layer_names), "."
  )
}

feature_universe <- rownames(post_obj)

extract_aligned_count_layer <- function(
    layer_name
) {
  
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
      "Unexpected features in layer: ",
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
  
  layer_matrix[
    feature_universe,
    ,
    drop = FALSE
  ]
}

cat(
  "Five count layers verified; alignment function defined.\n"
)

# ============================================================
# 7. Construct the aligned raw-count matrix
# ============================================================

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

stopifnot(
  identical(
    common_cells,
    colnames(reference_obj)
  )
)

cat(
  "Aligned raw-count matrix:",
  paste(dim(raw_counts), collapse = " x "),
  "\n"
)





# ============================================================
# 8. Define full-data HVG sensitivity function
# ============================================================

rm(count_matrix_list)
invisible(gc())

run_HVG_sensitivity_test <- function(
    hvg_count
) {
  
  test_start <- Sys.time()
  
  hvg_metadata <- metadata_source
  
  hvg_metadata$donor_id_analysis <- as.character(
    hvg_metadata[, "sample_id"]
  )
  
  test_obj <- CreateSeuratObject(
    counts = raw_counts,
    meta.data = hvg_metadata,
    assay = "RNA",
    project = paste0(
      "HVG_Sensitivity_",
      hvg_count
    ),
    min.cells = 0,
    min.features = 0
  )
  
  DefaultAssay(test_obj) <- "RNA"
  
  test_obj[["RNA"]] <- split(
    test_obj[["RNA"]],
    f = test_obj$donor_id_analysis
  )
  
  set.seed(42)
  
  test_obj <- NormalizeData(
    object = test_obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  
  test_obj <- FindVariableFeatures(
    object = test_obj,
    selection.method = "vst",
    nfeatures = hvg_count,
    verbose = FALSE
  )
  
  test_obj <- ScaleData(
    object = test_obj,
    features = VariableFeatures(test_obj),
    verbose = FALSE
  )
  
  test_obj <- RunPCA(
    object = test_obj,
    features = VariableFeatures(test_obj),
    npcs = baseline_PCs,
    seed.use = 42,
    verbose = FALSE
  )
  
  set.seed(42)
  
  test_obj <- IntegrateLayers(
    object = test_obj,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca",
    dims = seq_len(baseline_PCs),
    verbose = FALSE
  )
  
  test_obj <- FindNeighbors(
    object = test_obj,
    reduction = "integrated.rpca",
    dims = seq_len(baseline_PCs),
    k.param = baseline_k,
    verbose = FALSE
  )
  
  set.seed(clustering_seed)
  
  test_obj <- FindClusters(
    object = test_obj,
    resolution = fixed_resolution,
    algorithm = 3,
    random.seed = clustering_seed,
    verbose = FALSE
  )
  
  test_labels <- as.character(
    Idents(test_obj)
  )
  
  names(test_labels) <- colnames(test_obj)
  
  aligned_test_labels <- test_labels[
    names(reference_labels)
  ]
  
  if (anyNA(aligned_test_labels)) {
    stop(
      "Cell alignment failed for HVG count ",
      hvg_count
    )
  }
  
  test_runtime <- as.numeric(
    difftime(
      Sys.time(),
      test_start,
      units = "mins"
    )
  )
  
  result <- data.frame(
    parameter = "HVGs",
    value = hvg_count,
    HVGs = hvg_count,
    PCs = baseline_PCs,
    neighbour_k = baseline_k,
    clusters_recovered = length(
      unique(aligned_test_labels)
    ),
    ARI_vs_retained_solution =
      mclust::adjustedRandIndex(
        reference_labels,
        aligned_test_labels
      ),
    NMI_vs_retained_solution =
      aricode::NMI(
        reference_labels,
        aligned_test_labels
      ),
    runtime_seconds =
      test_runtime * 60,
    stringsAsFactors = FALSE
  )
  
  rm(test_obj)
  invisible(gc())
  
  result
}

cat("Full-data HVG sensitivity function defined.\n")


# ============================================================
# 9. Run 1000-HVG sensitivity test
# ============================================================

cat("Running full-data RPCA with 1000 HVGs...\n")

hvg_1000_result <-
  run_HVG_sensitivity_test(
    hvg_count = 1000
  )

print(
  hvg_1000_result %>%
    mutate(
      ARI_vs_retained_solution = round(
        ARI_vs_retained_solution,
        4
      ),
      NMI_vs_retained_solution = round(
        NMI_vs_retained_solution,
        4
      ),
      runtime_seconds = round(
        runtime_seconds,
        2
      )
    ),
  row.names = FALSE
)

write.csv(
  hvg_1000_result,
  file.path(
    output_dir,
    "checkpoint_HVG_1000_result.csv"
  ),
  row.names = FALSE
)

cat("1000-HVG test completed successfully.\n")


# ============================================================
# 10. Run 3000-HVG sensitivity test
# ============================================================

cat("Running full-data RPCA with 3000 HVGs...\n")

hvg_3000_result <-
  run_HVG_sensitivity_test(
    hvg_count = 3000
  )

print(
  hvg_3000_result %>%
    mutate(
      ARI_vs_retained_solution = round(
        ARI_vs_retained_solution,
        4
      ),
      NMI_vs_retained_solution = round(
        NMI_vs_retained_solution,
        4
      ),
      runtime_seconds = round(
        runtime_seconds,
        2
      )
    ),
  row.names = FALSE
)

write.csv(
  hvg_3000_result,
  file.path(
    output_dir,
    "checkpoint_HVG_3000_result.csv"
  ),
  row.names = FALSE
)

cat("3000-HVG test completed successfully.\n")

# ============================================================
# 11. Combine all parameter-sensitivity results
# ============================================================

hvg_baseline_result <- data.frame(
  parameter = "HVGs",
  value = baseline_HVGs,
  HVGs = baseline_HVGs,
  PCs = baseline_PCs,
  neighbour_k = baseline_k,
  clusters_recovered =
    reference_cluster_count,
  ARI_vs_retained_solution = 1,
  NMI_vs_retained_solution = 1,
  runtime_seconds = NA_real_,
  stringsAsFactors = FALSE
)

hvg_results <- bind_rows(
  hvg_1000_result,
  hvg_baseline_result,
  hvg_3000_result
)

all_parameter_results <- bind_rows(
  graph_sensitivity_results,
  hvg_results
) %>%
  mutate(
    parameter = factor(
      parameter,
      levels = c(
        "HVGs",
        "PCs",
        "Neighbour_k"
      )
    ),
    is_retained_setting = case_when(
      parameter == "HVGs" &
        value == baseline_HVGs ~ TRUE,
      parameter == "PCs" &
        value == baseline_PCs ~ TRUE,
      parameter == "Neighbour_k" &
        value == baseline_k ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  arrange(
    parameter,
    value
  )

stopifnot(
  nrow(all_parameter_results) == 9
)

write.csv(
  all_parameter_results,
  file.path(
    output_dir,
    "Table_all_parameter_sensitivity_results.csv"
  ),
  row.names = FALSE
)

print(
  all_parameter_results %>%
    mutate(
      ARI_vs_retained_solution = round(
        ARI_vs_retained_solution,
        4
      ),
      NMI_vs_retained_solution = round(
        NMI_vs_retained_solution,
        4
      ),
      runtime_seconds = round(
        runtime_seconds,
        2
      )
    ),
  row.names = FALSE
)




# ============================================================
# 12. Generate parameter-sensitivity figure
# ============================================================

agreement_long <- all_parameter_results %>%
  select(
    parameter,
    value,
    ARI_vs_retained_solution,
    NMI_vs_retained_solution
  ) %>%
  pivot_longer(
    cols = c(
      ARI_vs_retained_solution,
      NMI_vs_retained_solution
    ),
    names_to = "agreement_metric",
    values_to = "agreement"
  ) %>%
  mutate(
    agreement_metric = factor(
      agreement_metric,
      levels = c(
        "ARI_vs_retained_solution",
        "NMI_vs_retained_solution"
      ),
      labels = c(
        "ARI",
        "NMI"
      )
    )
  )

baseline_lines <- data.frame(
  parameter = factor(
    c(
      "HVGs",
      "PCs",
      "Neighbour_k"
    ),
    levels = levels(
      all_parameter_results$parameter
    )
  ),
  retained_value = c(
    baseline_HVGs,
    baseline_PCs,
    baseline_k
  )
)

parameter_sensitivity_figure <- ggplot(
  agreement_long,
  aes(
    x = value,
    y = agreement,
    colour = agreement_metric,
    group = agreement_metric
  )
) +
  geom_vline(
    data = baseline_lines,
    aes(xintercept = retained_value),
    inherit.aes = FALSE,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 3
  ) +
  facet_wrap(
    ~ parameter,
    scales = "free_x",
    nrow = 1,
    labeller = as_labeller(
      c(
        "HVGs" = "Highly variable genes",
        "PCs" = "Principal components",
        "Neighbour_k" = "Neighbour k"
      )
    )
  ) +
  scale_colour_manual(
    values = c(
      "ARI" = "#0072B2",
      "NMI" = "#D55E00"
    )
  ) +
  coord_cartesian(
    ylim = c(0.60, 1.02)
  ) +
  labs(
    title = "One-factor-at-a-time parameter sensitivity",
    subtitle = "Dashed lines indicate the retained parameter settings",
    x = "Parameter value",
    y = "Agreement with retained resolution-0.6 solution",
    colour = "Metric"
  ) +
  theme_bw(base_size = 13) +
  theme(
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
    legend.position = "bottom"
  )

print(parameter_sensitivity_figure)

ggsave(
  file.path(
    output_dir,
    "Figure_parameter_sensitivity_agreement.png"
  ),
  parameter_sensitivity_figure,
  width = 12,
  height = 4.8,
  dpi = 400,
  bg = "white"
)

ggsave(
  file.path(
    output_dir,
    "Figure_parameter_sensitivity_agreement.pdf"
  ),
  parameter_sensitivity_figure,
  width = 12,
  height = 4.8,
  bg = "white"
)

cat("Parameter-sensitivity figure saved.\n")



# ============================================================
# 13. Summarize non-baseline parameter robustness
# ============================================================

parameter_robustness_summary <-
  all_parameter_results %>%
  filter(!is_retained_setting) %>%
  group_by(parameter) %>%
  summarise(
    alternative_values = paste(
      value,
      collapse = ", "
    ),
    minimum_ARI = round(
      min(ARI_vs_retained_solution),
      4
    ),
    maximum_ARI = round(
      max(ARI_vs_retained_solution),
      4
    ),
    minimum_NMI = round(
      min(NMI_vs_retained_solution),
      4
    ),
    maximum_NMI = round(
      max(NMI_vs_retained_solution),
      4
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
  parameter_robustness_summary,
  file.path(
    output_dir,
    "Table_parameter_robustness_summary.csv"
  ),
  row.names = FALSE
)

print(
  as.data.frame(
    parameter_robustness_summary
  ),
  row.names = FALSE
)



# ============================================================
# 14. Save summary and session information
# ============================================================

summary_file <- file.path(
  output_dir,
  "parameter_sensitivity_summary.txt"
)

capture.output(
  {
    cat(
      "PARAMETER SENSITIVITY ANALYSIS\n",
      "==============================\n\n"
    )
    
    cat(
      "Retained configuration:\n",
      "- Highly variable genes: 2000\n",
      "- Principal components: 20\n",
      "- Neighbour k: 20\n",
      "- Clustering resolution: 0.6\n\n"
    )
    
    cat("ALL PARAMETER RESULTS\n")
    print(
      as.data.frame(all_parameter_results),
      row.names = FALSE
    )
    
    cat("\nALTERNATIVE-SETTING ROBUSTNESS\n")
    print(
      as.data.frame(
        parameter_robustness_summary
      ),
      row.names = FALSE
    )
    
    cat(
      "\nINTERPRETATION\n",
      "The retained parameter setting has an agreement ",
      "value of 1 by definition because it is the retained ",
      "reference solution; this is not treated as independent ",
      "evidence of superiority. Nearby alternative settings ",
      "retained moderate-to-high clustering agreement, with ",
      "ARI values from 0.6402 to 0.7837 and NMI values from ",
      "0.7514 to 0.8390. Cluster counts varied from 16 to 20. ",
      "The selected values of 2000 HVGs, 20 PCs, and k = 20 ",
      "represent intermediate settings that recovered 18 ",
      "clusters, whereas smaller or larger settings produced ",
      "moderate merging or fragmentation. Thus, broad cellular ",
      "structure remained reproducible, although fine-grained ",
      "cluster partitioning remained parameter-sensitive.\n",
      sep = ""
    )
  },
  file = summary_file
)

session_info_file <- file.path(
  output_dir,
  "parameter_sensitivity_session_info.txt"
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