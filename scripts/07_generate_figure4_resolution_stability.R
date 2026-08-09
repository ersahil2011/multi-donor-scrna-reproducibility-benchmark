# ============================================================
# 07_generate_figure4_resolution_stability.R
#
# Purpose:
# Generate Figure 4:
# (a) Number of clusters across resolutions
# (b) ARI and NMI agreement with reference resolution 0.6
# ============================================================

required_packages <- c("ggplot2", "tidyr", "dplyr", "patchwork")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(patchwork)
})

# ------------------------------------------------------------
# 1. Resolution stability data
# ------------------------------------------------------------

stability_df <- data.frame(
  resolution = c(0.2, 0.4, 0.5, 0.6, 0.8, 1.0),
  clusters = c(10, 14, 18, 18, 19, 21),
  ARI_vs_0_6 = c(0.4866, 0.6620, 0.6676, 1.0000, 0.8800, 0.8210),
  NMI_vs_0_6 = c(0.6581, 0.7829, 0.7731, 1.0000, 0.9045, 0.8631)
)

output_dir <- "figure4_resolution_stability_outputs"

if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}

write.csv(
  stability_df,
  file.path(output_dir, "Figure4_resolution_stability_values.csv"),
  row.names = FALSE
)

# ------------------------------------------------------------
# 2. Panel A: Cluster count across resolutions
# ------------------------------------------------------------

p_clusters <- ggplot(
  stability_df,
  aes(x = resolution, y = clusters)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  scale_x_continuous(
    breaks = stability_df$resolution
  ) +
  scale_y_continuous(
    breaks = seq(10, 22, by = 2)
  ) +
  labs(
    title = "(a) Cluster count across resolutions",
    x = "Clustering resolution",
    y = "Number of clusters"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0)
  )

# ------------------------------------------------------------
# 3. Panel B: ARI and NMI agreement
# ------------------------------------------------------------

agreement_df <- stability_df %>%
  select(resolution, ARI_vs_0_6, NMI_vs_0_6) %>%
  pivot_longer(
    cols = c(ARI_vs_0_6, NMI_vs_0_6),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    metric = recode(
      metric,
      ARI_vs_0_6 = "ARI vs resolution 0.6",
      NMI_vs_0_6 = "NMI vs resolution 0.6"
    )
  )

p_agreement <- ggplot(
  agreement_df,
  aes(x = resolution, y = value, shape = metric, linetype = metric)
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.8) +
  scale_x_continuous(
    breaks = stability_df$resolution
  ) +
  scale_y_continuous(
    limits = c(0, 1.05),
    breaks = seq(0, 1, by = 0.2)
  ) +
  labs(
    title = "(b) Agreement with reference solution",
    x = "Clustering resolution",
    y = "Agreement score",
    shape = "Metric",
    linetype = "Metric"
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0),
    legend.title = element_text(face = "bold"),
    legend.position = "right"
  )

# ------------------------------------------------------------
# 4. Combine panels
# ------------------------------------------------------------

figure4 <- p_clusters + p_agreement + plot_layout(ncol = 2)

# ------------------------------------------------------------
# 5. Save outputs
# ------------------------------------------------------------

ggsave(
  filename = file.path(output_dir, "Figure4_resolution_wise_clustering_stability.png"),
  plot = figure4,
  width = 12,
  height = 5,
  dpi = 600
)

ggsave(
  filename = file.path(output_dir, "Figure4_resolution_wise_clustering_stability.pdf"),
  plot = figure4,
  width = 12,
  height = 5
)

cat("\nFigure 4 generated successfully.\n")
cat("Output folder:", output_dir, "\n")
cat("Files created:\n")
cat("1. Figure4_resolution_wise_clustering_stability.png\n")
cat("2. Figure4_resolution_wise_clustering_stability.pdf\n")