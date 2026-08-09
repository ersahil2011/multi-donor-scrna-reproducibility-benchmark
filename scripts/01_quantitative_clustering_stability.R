# ============================================================
# 01_quantitative_clustering_stability.R
#
# Purpose:
# Quantitatively evaluate clustering stability across multiple
# graph-clustering resolution settings using the final RPCA-
# integrated Seurat object.
#
# Outputs:
# 1. Resolution-wise cluster counts
# 2. ARI and NMI versus reference resolution 0.6
# 3. Silhouette summary on a 5,000-cell reproducible subset
# 4. Stability plots and session information
# ============================================================


# ============================================================
# 1. Package setup
# ============================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "dplyr",
  "ggplot2",
  "tidyr",
  "mclust",
  "aricode",
  "cluster",
  "Matrix",
  "patchwork"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(mclust)
  library(aricode)
  library(cluster)
  library(Matrix)
  library(patchwork)
})


# ============================================================
# 2. User parameters
# ============================================================

input_file <- "sce_final_clusters_locked_rpca_pcs20_res06.rds"

output_dir <- "stability_outputs"

reference_resolution <- 0.6

evaluated_resolutions <- c(0.2, 0.4, 0.5, 0.6, 0.8, 1.0)

dims_use <- 1:20

silhouette_sample_size <- 5000

silhouette_seed <- 42


# ============================================================
# 3. Create output directory
# ============================================================

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}


# ============================================================
# 4. Load object
# ============================================================

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

cat("\nLoading object:\n")
cat(input_file, "\n")

obj <- readRDS(input_file)

cat("\nObject dimensions:\n")
print(dim(obj))

cat("\nAvailable assays:\n")
print(Assays(obj))

cat("\nAvailable reductions:\n")
print(Reductions(obj))

cat("\nAvailable metadata columns:\n")
print(colnames(obj@meta.data))


# ============================================================
# 5. Helper functions
# ============================================================

format_resolution <- function(x) {
  # Converts 0.6 to "0.6" and also supports Seurat metadata pattern matching
  as.character(x)
}

find_cluster_column <- function(metadata_cols, resolution_value) {
  
  resolution_text <- format_resolution(resolution_value)
  
  candidate_patterns <- c(
    paste0("res\\.", resolution_text, "$"),
    paste0("res_", gsub("\\.", "_", resolution_text), "$"),
    paste0("res", gsub("\\.", "", resolution_text), "$"),
    paste0("snn_res\\.", resolution_text, "$"),
    paste0("RNA_snn_res\\.", resolution_text, "$"),
    paste0("integrated_snn_res\\.", resolution_text, "$"),
    paste0("rpca_snn_res\\.", resolution_text, "$"),
    paste0("SCT_snn_res\\.", resolution_text, "$")
  )
  
  for (pattern in candidate_patterns) {
    matched <- grep(
      pattern,
      metadata_cols,
      value = TRUE,
      ignore.case = TRUE
    )
    
    if (length(matched) > 0) {
      return(matched[1])
    }
  }
  
  # Broader fallback
  broad_pattern <- paste0("res.*", resolution_text)
  
  matched <- grep(
    broad_pattern,
    metadata_cols,
    value = TRUE,
    ignore.case = TRUE
  )
  
  if (length(matched) > 0) {
    return(matched[1])
  }
  
  return(NA_character_)
}

detect_reduction_for_embedding <- function(obj) {
  
  reductions_available <- Reductions(obj)
  
  preferred <- c(
    "pca",
    "integrated_pca",
    "rpca",
    "harmony"
  )
  
  for (red in preferred) {
    if (red %in% reductions_available) {
      return(red)
    }
  }
  
  pca_like <- reductions_available[
    grepl("pca|rpca|integrated", reductions_available, ignore.case = TRUE)
  ]
  
  if (length(pca_like) > 0) {
    return(pca_like[1])
  }
  
  stop(
    "Could not detect PCA/RPCA-like reduction for silhouette analysis. Available reductions: ",
    paste(reductions_available, collapse = ", ")
  )
}


# ============================================================
# 6. Detect cluster columns for all evaluated resolutions
# ============================================================

metadata_cols <- colnames(obj@meta.data)

cluster_column_map <- data.frame(
  resolution = evaluated_resolutions,
  cluster_column = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_along(evaluated_resolutions)) {
  
  res <- evaluated_resolutions[i]
  
  detected_col <- find_cluster_column(
    metadata_cols = metadata_cols,
    resolution_value = res
  )
  
  cluster_column_map$cluster_column[i] <- detected_col
}

cat("\nDetected cluster columns:\n")
print(cluster_column_map)

if (any(is.na(cluster_column_map$cluster_column))) {
  missing_res <- cluster_column_map$resolution[
    is.na(cluster_column_map$cluster_column)
  ]
  
  stop(
    "Could not detect cluster columns for resolution(s): ",
    paste(missing_res, collapse = ", "),
    "\nAvailable metadata columns are:\n",
    paste(metadata_cols, collapse = ", ")
  )
}

write.csv(
  cluster_column_map,
  file.path(output_dir, "Table_detected_cluster_columns.csv"),
  row.names = FALSE
)


# ============================================================
# 7. Extract reference clustering at resolution 0.6
# ============================================================

reference_col <- cluster_column_map$cluster_column[
  cluster_column_map$resolution == reference_resolution
]

if (length(reference_col) != 1 || is.na(reference_col)) {
  stop("Reference resolution column not detected properly.")
}

reference_labels <- as.factor(obj@meta.data[[reference_col]])

cat("\nReference column used:\n")
cat(reference_col, "\n")

cat("\nReference cluster count:\n")
print(length(unique(reference_labels)))


# ============================================================
# 8. Compute ARI/NMI and cluster counts across resolutions
# ============================================================

stability_results <- lapply(
  seq_len(nrow(cluster_column_map)),
  function(i) {
    
    res <- cluster_column_map$resolution[i]
    col <- cluster_column_map$cluster_column[i]
    
    current_labels <- as.factor(obj@meta.data[[col]])
    
    data.frame(
      resolution = res,
      cluster_column = col,
      number_of_clusters = length(unique(current_labels)),
      ARI_vs_resolution_0_6 = mclust::adjustedRandIndex(
        reference_labels,
        current_labels
      ),
      NMI_vs_resolution_0_6 = aricode::NMI(
        reference_labels,
        current_labels
      ),
      stringsAsFactors = FALSE
    )
  }
) %>%
  bind_rows()

cat("\nResolution-wise stability results:\n")
print(stability_results)

write.csv(
  stability_results,
  file.path(output_dir, "Table_resolution_wise_clustering_stability.csv"),
  row.names = FALSE
)


# ============================================================
# 9. Generate stability plots
# ============================================================

p_cluster_count <- ggplot(
  stability_results,
  aes(
    x = resolution,
    y = number_of_clusters
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  scale_x_continuous(
    breaks = evaluated_resolutions
  ) +
  scale_y_continuous(
    breaks = seq(
      min(stability_results$number_of_clusters),
      max(stability_results$number_of_clusters),
      by = 1
    )
  ) +
  labs(
    title = "Cluster count across resolutions",
    x = "Clustering resolution",
    y = "Number of clusters"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0)
  )

agreement_long <- stability_results %>%
  select(
    resolution,
    ARI_vs_resolution_0_6,
    NMI_vs_resolution_0_6
  ) %>%
  pivot_longer(
    cols = c(
      ARI_vs_resolution_0_6,
      NMI_vs_resolution_0_6
    ),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      ARI_vs_resolution_0_6 = "ARI",
      NMI_vs_resolution_0_6 = "NMI"
    )
  )

p_agreement <- ggplot(
  agreement_long,
  aes(
    x = resolution,
    y = value,
    shape = metric,
    linetype = metric
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  scale_x_continuous(
    breaks = evaluated_resolutions
  ) +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(
    title = "Agreement with resolution 0.6 reference",
    x = "Clustering resolution",
    y = "Agreement score",
    shape = "Metric",
    linetype = "Metric"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    legend.title = element_text(face = "bold")
  )

figure_stability <- p_cluster_count + p_agreement + patchwork::plot_layout(ncol = 2)

ggsave(
  filename = file.path(output_dir, "Figure_resolution_stability_cluster_count_and_agreement.png"),
  plot = figure_stability,
  width = 12,
  height = 5,
  dpi = 600,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Figure_resolution_stability_cluster_count_and_agreement.pdf"),
  plot = figure_stability,
  width = 12,
  height = 5,
  limitsize = FALSE
)


# ============================================================
# 10. Silhouette analysis on reproducible subset
# ============================================================

cat("\nStarting silhouette analysis...\n")

embedding_reduction <- detect_reduction_for_embedding(obj)

cat("\nEmbedding reduction used for silhouette:\n")
cat(embedding_reduction, "\n")

embedding_matrix <- Embeddings(obj, reduction = embedding_reduction)

if (ncol(embedding_matrix) < max(dims_use)) {
  stop(
    "Embedding has fewer dimensions than required. Available dimensions: ",
    ncol(embedding_matrix)
  )
}

embedding_matrix <- embedding_matrix[, dims_use, drop = FALSE]

set.seed(silhouette_seed)

if (nrow(embedding_matrix) > silhouette_sample_size) {
  sampled_cells <- sample(
    rownames(embedding_matrix),
    size = silhouette_sample_size,
    replace = FALSE
  )
} else {
  sampled_cells <- rownames(embedding_matrix)
}

embedding_sample <- embedding_matrix[sampled_cells, , drop = FALSE]

cluster_sample <- as.factor(reference_labels[sampled_cells])

cat("\nSilhouette sampled cells:\n")
cat(length(sampled_cells), "\n")

cat("\nComputing distance matrix for silhouette subset...\n")

distance_matrix <- dist(embedding_sample)

silhouette_object <- cluster::silhouette(
  as.integer(cluster_sample),
  distance_matrix
)

silhouette_df <- as.data.frame(silhouette_object)

silhouette_summary <- data.frame(
  reference_resolution = reference_resolution,
  reference_cluster_column = reference_col,
  embedding_reduction = embedding_reduction,
  dimensions_used = paste(dims_use, collapse = ","),
  sampled_cells = length(sampled_cells),
  seed = silhouette_seed,
  mean_silhouette_width = mean(silhouette_df$sil_width),
  median_silhouette_width = median(silhouette_df$sil_width),
  minimum_silhouette_width = min(silhouette_df$sil_width),
  maximum_silhouette_width = max(silhouette_df$sil_width),
  stringsAsFactors = FALSE
)

cat("\nSilhouette summary:\n")
print(silhouette_summary)

write.csv(
  silhouette_summary,
  file.path(output_dir, "Table_silhouette_summary_resolution_0_6.csv"),
  row.names = FALSE
)

write.csv(
  data.frame(
    cell_barcode = sampled_cells,
    reference_cluster = cluster_sample,
    silhouette_width = silhouette_df$sil_width,
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "Table_silhouette_cell_level_values.csv"),
  row.names = FALSE
)

p_silhouette <- ggplot(
  silhouette_df,
  aes(x = sil_width)
) +
  geom_histogram(
    bins = 40
  ) +
  labs(
    title = "Silhouette width distribution at resolution 0.6",
    x = "Silhouette width",
    y = "Number of sampled cells"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0)
  )

ggsave(
  filename = file.path(output_dir, "Figure_silhouette_width_distribution.png"),
  plot = p_silhouette,
  width = 7,
  height = 5,
  dpi = 600,
  limitsize = FALSE
)

ggsave(
  filename = file.path(output_dir, "Figure_silhouette_width_distribution.pdf"),
  plot = p_silhouette,
  width = 7,
  height = 5,
  limitsize = FALSE
)


# ============================================================
# 11. Save text summary
# ============================================================

summary_text <- c(
  "Quantitative clustering stability analysis",
  "==========================================",
  paste0("Input object: ", input_file),
  paste0("Reference resolution: ", reference_resolution),
  paste0("Reference cluster column: ", reference_col),
  paste0("Evaluated resolutions: ", paste(evaluated_resolutions, collapse = ", ")),
  paste0("Embedding reduction used for silhouette: ", embedding_reduction),
  paste0("Dimensions used for silhouette: ", paste(dims_use, collapse = ",")),
  paste0("Silhouette sample size: ", length(sampled_cells)),
  paste0("Silhouette seed: ", silhouette_seed),
  paste0("Mean silhouette width: ", round(silhouette_summary$mean_silhouette_width, 4)),
  paste0("Median silhouette width: ", round(silhouette_summary$median_silhouette_width, 4))
)

writeLines(
  summary_text,
  file.path(output_dir, "clustering_stability_summary.txt")
)


# ============================================================
# 12. Save session information
# ============================================================

sink(file.path(output_dir, "clustering_stability_session_info.txt"))
print(sessionInfo())
sink()


# ============================================================
# 13. Final message
# ============================================================

cat("\n====================================================\n")
cat("QUANTITATIVE CLUSTERING STABILITY ANALYSIS COMPLETED\n")
cat("====================================================\n")
cat("Output folder:", output_dir, "\n\n")

cat("Files created:\n")
cat("1. Table_detected_cluster_columns.csv\n")
cat("2. Table_resolution_wise_clustering_stability.csv\n")
cat("3. Table_silhouette_summary_resolution_0_6.csv\n")
cat("4. Table_silhouette_cell_level_values.csv\n")
cat("5. Figure_resolution_stability_cluster_count_and_agreement.png\n")
cat("6. Figure_resolution_stability_cluster_count_and_agreement.pdf\n")
cat("7. Figure_silhouette_width_distribution.png\n")
cat("8. Figure_silhouette_width_distribution.pdf\n")
cat("9. clustering_stability_summary.txt\n")
cat("10. clustering_stability_session_info.txt\n\n")

cat("Resolution-wise stability table:\n")
print(stability_results)

cat("\nSilhouette summary:\n")
print(silhouette_summary)