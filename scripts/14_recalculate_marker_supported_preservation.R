# ============================================================
# 1. Verify required saved inputs
# ============================================================

rm(list = ls())
gc()

input_files <- c(
  cell_cluster_map =
    "cell_labels_final.csv",
  
  balanced_cell_results =
    file.path(
      "balanced_celltype_preservation_outputs",
      "cell_level_label_preservation_all_iterations.csv"
    ),
  
  marker_annotation_map =
    file.path(
      "reviewer_revision_outputs",
      "10_marker_annotation_validation",
      "TableS2_marker_supported_annotation_review.csv"
    )
)

print(file.exists(input_files))



# ============================================================
# 2. Load packages and saved results
# ============================================================

required_packages <- c(
  "dplyr",
  "tidyr",
  "mclust",
  "aricode",
  "ggplot2"
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
  library(dplyr)
  library(tidyr)
  library(mclust)
  library(aricode)
  library(ggplot2)
})

output_dir <- file.path(
  "reviewer_revision_outputs",
  "14_marker_supported_preservation"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cell_cluster_map <- read.csv(
  unname(input_files["cell_cluster_map"]),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

balanced_cells <- read.csv(
  unname(input_files["balanced_cell_results"]),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

marker_annotation_map <- read.csv(
  unname(input_files["marker_annotation_map"]),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

stopifnot(
  all(c("cell", "cluster") %in% colnames(cell_cluster_map)),
  all(
    c(
      "iteration",
      "seed",
      "cell_barcode",
      "donor_id",
      "balanced_cluster"
    ) %in% colnames(balanced_cells)
  ),
  all(
    c(
      "cluster",
      "marker_supported_label"
    ) %in% colnames(marker_annotation_map)
  )
)

cat(
  "Cell-cluster map:",
  nrow(cell_cluster_map),
  "rows.\n"
)

cat(
  "Saved balanced cell records:",
  nrow(balanced_cells),
  "rows.\n"
)

cat(
  "Marker-supported cluster annotations:",
  nrow(marker_annotation_map),
  "rows.\n"
)

cat("Output folder:", output_dir, "\n")


# ============================================================
# 3. Attach marker-supported labels to sampled cells
# ============================================================

cell_reference_new <- cell_cluster_map %>%
  transmute(
    cell_barcode = as.character(cell),
    full_data_cluster = as.character(cluster)
  ) %>%
  left_join(
    marker_annotation_map %>%
      transmute(
        full_data_cluster = as.character(cluster),
        marker_supported_reference =
          as.character(marker_supported_label)
      ),
    by = "full_data_cluster",
    relationship = "many-to-one"
  )

stopifnot(
  nrow(cell_reference_new) == 43110,
  !anyDuplicated(cell_reference_new$cell_barcode),
  !anyNA(cell_reference_new$marker_supported_reference)
)

balanced_marker_cells <- balanced_cells %>%
  transmute(
    iteration = as.integer(iteration),
    seed = as.integer(seed),
    cell_barcode = as.character(cell_barcode),
    donor_id = as.character(donor_id),
    balanced_cluster = as.character(balanced_cluster)
  ) %>%
  left_join(
    cell_reference_new,
    by = "cell_barcode",
    relationship = "many-to-one"
  )

iteration_sizes <- balanced_marker_cells %>%
  count(
    iteration,
    name = "sampled_cells"
  )

stopifnot(
  nrow(balanced_marker_cells) == 347750,
  !anyNA(balanced_marker_cells$full_data_cluster),
  !anyNA(
    balanced_marker_cells$marker_supported_reference
  ),
  nrow(iteration_sizes) == 50,
  all(iteration_sizes$sampled_cells == 6955)
)

cat(
  "Marker-supported labels attached to",
  nrow(balanced_marker_cells),
  "saved cell records.\n"
)

cat(
  "Verified:",
  nrow(iteration_sizes),
  "iterations with",
  unique(iteration_sizes$sampled_cells),
  "cells each.\n"
)

print(
  table(
    cell_reference_new$marker_supported_reference
  )
)





# ============================================================
# 4. Assign majority marker-supported label to each rerun cluster
# ============================================================

cluster_label_counts <- balanced_marker_cells %>%
  count(
    iteration,
    seed,
    balanced_cluster,
    marker_supported_reference,
    name = "cells_with_label"
  ) %>%
  group_by(
    iteration,
    seed,
    balanced_cluster
  ) %>%
  mutate(
    cluster_size = sum(cells_with_label),
    cluster_purity = cells_with_label / cluster_size
  ) %>%
  ungroup()

marker_supported_cluster_mapping <-
  cluster_label_counts %>%
  arrange(
    iteration,
    balanced_cluster,
    desc(cells_with_label),
    marker_supported_reference
  ) %>%
  group_by(
    iteration,
    seed,
    balanced_cluster
  ) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  transmute(
    iteration,
    seed,
    balanced_cluster,
    predicted_marker_supported_label =
      marker_supported_reference,
    cluster_size,
    cluster_purity
  )

balanced_marker_predictions <-
  balanced_marker_cells %>%
  left_join(
    marker_supported_cluster_mapping,
    by = c(
      "iteration",
      "seed",
      "balanced_cluster"
    ),
    relationship = "many-to-one"
  ) %>%
  mutate(
    marker_label_preserved =
      marker_supported_reference ==
      predicted_marker_supported_label
  )

stopifnot(
  nrow(balanced_marker_predictions) == 347750,
  !anyNA(
    balanced_marker_predictions$
      predicted_marker_supported_label
  ),
  !anyNA(
    balanced_marker_predictions$
      marker_label_preserved
  )
)

mapping_check <- marker_supported_cluster_mapping %>%
  group_by(iteration) %>%
  summarise(
    clusters_mapped = n(),
    mean_cluster_purity = mean(cluster_purity),
    minimum_cluster_purity = min(cluster_purity),
    .groups = "drop"
  )

cat(
  "Majority labels assigned to",
  nrow(marker_supported_cluster_mapping),
  "iteration-specific clusters.\n"
)

cat(
  "Clusters mapped per iteration:",
  min(mapping_check$clusters_mapped),
  "to",
  max(mapping_check$clusters_mapped),
  "\n"
)

cat(
  "No missing predicted labels:",
  !anyNA(
    balanced_marker_predictions$
      predicted_marker_supported_label
  ),
  "\n"
)


# ============================================================
# 5. Calculate marker-supported preservation metrics
# ============================================================

all_marker_labels <- sort(
  unique(
    cell_reference_new$marker_supported_reference
  )
)

calculate_macro_F1 <- function(
    reference_labels,
    predicted_labels,
    label_levels
) {
  
  F1_values <- vapply(
    label_levels,
    function(current_label) {
      
      true_positive <- sum(
        reference_labels == current_label &
          predicted_labels == current_label
      )
      
      false_positive <- sum(
        reference_labels != current_label &
          predicted_labels == current_label
      )
      
      false_negative <- sum(
        reference_labels == current_label &
          predicted_labels != current_label
      )
      
      precision_denominator <-
        true_positive + false_positive
      
      recall_denominator <-
        true_positive + false_negative
      
      precision <- if (
        precision_denominator == 0
      ) {
        0
      } else {
        true_positive / precision_denominator
      }
      
      recall <- if (
        recall_denominator == 0
      ) {
        0
      } else {
        true_positive / recall_denominator
      }
      
      if ((precision + recall) == 0) {
        0
      } else {
        2 * precision * recall /
          (precision + recall)
      }
    },
    numeric(1)
  )
  
  mean(F1_values)
}

marker_supported_metrics <-
  balanced_marker_predictions %>%
  group_by(iteration, seed) %>%
  group_modify(
    function(iteration_data, iteration_key) {
      
      reference_labels <- as.character(
        iteration_data$marker_supported_reference
      )
      
      predicted_labels <- as.character(
        iteration_data$
          predicted_marker_supported_label
      )
      
      reference_codes <- as.integer(
        factor(
          reference_labels,
          levels = all_marker_labels
        )
      )
      
      predicted_codes <- as.integer(
        factor(
          predicted_labels,
          levels = all_marker_labels
        )
      )
      
      data.frame(
        cells_per_donor = min(
          table(iteration_data$donor_id)
        ),
        total_balanced_cells =
          nrow(iteration_data),
        clusters_recovered =
          n_distinct(
            iteration_data$balanced_cluster
          ),
        reference_populations_present =
          n_distinct(reference_labels),
        predicted_populations_recovered =
          n_distinct(predicted_labels),
        marker_label_agreement =
          mean(reference_labels == predicted_labels),
        macro_F1 = calculate_macro_F1(
          reference_labels,
          predicted_labels,
          all_marker_labels
        ),
        ARI_marker_structure =
          mclust::adjustedRandIndex(
            reference_codes,
            predicted_codes
          ),
        NMI_marker_structure =
          aricode::NMI(
            reference_codes,
            predicted_codes
          )
      )
    }
  ) %>%
  ungroup()

stopifnot(
  nrow(marker_supported_metrics) == 50,
  all(
    marker_supported_metrics$
      total_balanced_cells == 6955
  ),
  all(
    marker_supported_metrics$
      reference_populations_present == 5
  )
)

metric_summary <- marker_supported_metrics %>%
  summarise(
    mean_agreement = mean(marker_label_agreement),
    sd_agreement = sd(marker_label_agreement),
    mean_macro_F1 = mean(macro_F1),
    sd_macro_F1 = sd(macro_F1),
    mean_ARI = mean(ARI_marker_structure),
    sd_ARI = sd(ARI_marker_structure),
    mean_NMI = mean(NMI_marker_structure),
    sd_NMI = sd(NMI_marker_structure),
    minimum_populations_recovered =
      min(predicted_populations_recovered),
    maximum_populations_recovered =
      max(predicted_populations_recovered)
  )

print(
  metric_summary %>%
    mutate(
      across(
        where(is.numeric),
        function(x) round(x, 4)
      )
    ),
  row.names = FALSE
)


print(metric_summary, width = Inf)




# ============================================================
# 6. Calculate population-specific preservation
# ============================================================

population_metrics <-
  balanced_marker_predictions %>%
  group_by(iteration, seed) %>%
  group_modify(
    function(iteration_data, iteration_key) {
      
      bind_rows(
        lapply(
          all_marker_labels,
          function(current_label) {
            
            reference_labels <- as.character(
              iteration_data$
                marker_supported_reference
            )
            
            predicted_labels <- as.character(
              iteration_data$
                predicted_marker_supported_label
            )
            
            true_positive <- sum(
              reference_labels == current_label &
                predicted_labels == current_label
            )
            
            false_positive <- sum(
              reference_labels != current_label &
                predicted_labels == current_label
            )
            
            false_negative <- sum(
              reference_labels == current_label &
                predicted_labels != current_label
            )
            
            precision <- if (
              true_positive + false_positive == 0
            ) {
              0
            } else {
              true_positive /
                (true_positive + false_positive)
            }
            
            recall <- if (
              true_positive + false_negative == 0
            ) {
              0
            } else {
              true_positive /
                (true_positive + false_negative)
            }
            
            F1 <- if (
              precision + recall == 0
            ) {
              0
            } else {
              2 * precision * recall /
                (precision + recall)
            }
            
            data.frame(
              marker_supported_population =
                current_label,
              reference_cells = sum(
                reference_labels == current_label
              ),
              predicted_cells = sum(
                predicted_labels == current_label
              ),
              true_positive = true_positive,
              precision = precision,
              recall = recall,
              F1 = F1
            )
          }
        )
      )
    }
  ) %>%
  ungroup()

population_summary <- population_metrics %>%
  group_by(marker_supported_population) %>%
  summarise(
    mean_precision = mean(precision),
    sd_precision = sd(precision),
    mean_recall = mean(recall),
    sd_recall = sd(recall),
    mean_F1 = mean(F1),
    sd_F1 = sd(F1),
    minimum_F1 = min(F1),
    maximum_F1 = max(F1),
    .groups = "drop"
  )

stopifnot(
  nrow(population_metrics) == 250,
  nrow(population_summary) == 5
)

print(
  population_summary %>%
    mutate(
      across(
        where(is.numeric),
        function(x) round(x, 4)
      )
    ),
  width = Inf,
  row.names = FALSE
)


# ============================================================
# 7. Save corrected preservation tables
# ============================================================

full_data_population_distribution <-
  cell_reference_new %>%
  count(
    marker_supported_reference,
    name = "number_of_cells"
  ) %>%
  mutate(
    percentage_of_total_cells = round(
      100 * number_of_cells /
        sum(number_of_cells),
      2
    )
  ) %>%
  arrange(
    desc(number_of_cells)
  )

write.csv(
  marker_supported_metrics,
  file.path(
    output_dir,
    "Table_marker_supported_preservation_all_iterations.csv"
  ),
  row.names = FALSE
)

write.csv(
  metric_summary,
  file.path(
    output_dir,
    "Table_marker_supported_preservation_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  population_metrics,
  file.path(
    output_dir,
    "Table_marker_supported_population_metrics_all_iterations.csv"
  ),
  row.names = FALSE
)

write.csv(
  population_summary,
  file.path(
    output_dir,
    "Table_marker_supported_population_preservation_summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  marker_supported_cluster_mapping,
  file.path(
    output_dir,
    "marker_supported_cluster_mapping_all_iterations.csv"
  ),
  row.names = FALSE
)

write.csv(
  balanced_marker_predictions,
  file.path(
    output_dir,
    "cell_level_marker_supported_preservation_all_iterations.csv"
  ),
  row.names = FALSE
)

write.csv(
  full_data_population_distribution,
  file.path(
    output_dir,
    "full_data_marker_supported_population_distribution.csv"
  ),
  row.names = FALSE
)

required_output_files <- c(
  "Table_marker_supported_preservation_all_iterations.csv",
  "Table_marker_supported_preservation_summary.csv",
  "Table_marker_supported_population_metrics_all_iterations.csv",
  "Table_marker_supported_population_preservation_summary.csv",
  "marker_supported_cluster_mapping_all_iterations.csv",
  "cell_level_marker_supported_preservation_all_iterations.csv",
  "full_data_marker_supported_population_distribution.csv"
)

output_status <- file.exists(
  file.path(
    output_dir,
    required_output_files
  )
)

cat(
  "Corrected output tables saved:",
  sum(output_status),
  "of",
  length(output_status),
  "\n"
)

stopifnot(all(output_status))





# ============================================================
# 8. Generate marker-supported preservation figure
# ============================================================

overall_plot_data <- marker_supported_metrics %>%
  select(
    iteration,
    marker_label_agreement,
    macro_F1,
    ARI_marker_structure,
    NMI_marker_structure
  ) %>%
  pivot_longer(
    cols = -iteration,
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    panel = "A. Overall preservation",
    metric = recode(
      metric,
      marker_label_agreement = "Label agreement",
      macro_F1 = "Macro-F1",
      ARI_marker_structure = "ARI",
      NMI_marker_structure = "NMI"
    )
  )

population_plot_data <- population_metrics %>%
  transmute(
    iteration,
    metric = marker_supported_population,
    value = F1,
    panel = "B. Population-specific F1"
  )

metric_order <- c(
  "Label agreement",
  "Macro-F1",
  "ARI",
  "NMI",
  "Endothelial cells",
  "Fibroblast/stromal cells",
  "Mural/perivascular cells",
  "Myeloid cells",
  "T cells"
)

panel_order <- c(
  "A. Overall preservation",
  "B. Population-specific F1"
)

preservation_plot_data <- bind_rows(
  overall_plot_data,
  population_plot_data
) %>%
  mutate(
    metric = factor(
      metric,
      levels = rev(metric_order)
    ),
    panel = factor(
      panel,
      levels = panel_order
    )
  )

write.csv(
  preservation_plot_data,
  file.path(
    output_dir,
    "Figure_marker_supported_preservation_plot_data.csv"
  ),
  row.names = FALSE
)

preservation_figure <- ggplot(
  preservation_plot_data,
  aes(
    x = value,
    y = metric,
    fill = panel
  )
) +
  geom_boxplot(
    width = 0.58,
    outlier.shape = NA,
    alpha = 0.85,
    colour = "black",
    linewidth = 0.6,
    orientation = "y"
  ) +
  geom_jitter(
    width = 0,
    height = 0.10,
    size = 1.4,
    alpha = 0.55,
    colour = "black"
  ) +
  facet_wrap(
    facets = vars(panel),
    ncol = 1,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = c(
      "A. Overall preservation" = "#0072B2",
      "B. Population-specific F1" = "#D55E00"
    )
  ) +
  scale_x_continuous(
    limits = c(0, 1.02),
    breaks = seq(0, 1, 0.2),
    expand = expansion(
      mult = c(0, 0.02)
    )
  ) +
  labs(
    title =
      "Marker-supported preservation across reruns",
    subtitle =
      "Fifty donor-balanced RPCA reruns; 6,955 cells per iteration",
    x = "Metric value",
    y = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 12
    ),
    strip.background = element_rect(
      fill = "grey85",
      colour = "black"
    ),
    strip.text = element_text(
      face = "bold",
      size = 12
    ),
    axis.text.x = element_text(
      size = 10,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 10,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 12
    ),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.margin = margin(
      10,
      15,
      10,
      15
    )
  )

print(preservation_figure)

figure_png_file <- file.path(
  output_dir,
  "Figure_marker_supported_preservation.png"
)

figure_pdf_file <- file.path(
  output_dir,
  "Figure_marker_supported_preservation.pdf"
)

ggsave(
  filename = figure_png_file,
  plot = preservation_figure,
  width = 10,
  height = 8,
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = figure_pdf_file,
  plot = preservation_figure,
  width = 10,
  height = 8,
  bg = "white"
)

cat(
  "Corrected preservation figures saved:",
  file.exists(figure_png_file) &&
    file.exists(figure_pdf_file),
  "\n"
)




# ============================================================
# 9. Save summary, checksums, session information and manifest
# ============================================================

input_checksums <- data.frame(
  input_name = names(input_files),
  input_path = unname(input_files),
  MD5_checksum = unname(
    tools::md5sum(
      unname(input_files)
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  input_checksums,
  file.path(
    output_dir,
    "marker_supported_preservation_input_checksums.csv"
  ),
  row.names = FALSE
)

summary_file <- file.path(
  output_dir,
  "marker_supported_preservation_summary.txt"
)

capture.output(
  {
    cat(
      "MARKER-SUPPORTED CELL-POPULATION PRESERVATION\n"
    )
    cat(
      "================================================\n\n"
    )
    
    cat(
      "Method:\n",
      "Existing cell-level results from 50 donor-balanced ",
      "RPCA reruns were re-evaluated using the ",
      "marker-supported full-data annotations.\n",
      "No RPCA rerun or reclustering was performed during ",
      "this recalculation.\n\n",
      sep = ""
    )
    
    cat(
      "Iterations: 50\n",
      "Cells per iteration: 6,955\n",
      "Cells per donor: 1,391\n",
      "Marker-supported populations: 5\n\n",
      sep = ""
    )
    
    cat("Overall preservation summary:\n")
    print(
      as.data.frame(metric_summary),
      row.names = FALSE
    )
    
    cat("\nPopulation-specific preservation summary:\n")
    print(
      as.data.frame(population_summary),
      row.names = FALSE
    )
    
    cat("\nFull-data population distribution:\n")
    print(
      as.data.frame(
        full_data_population_distribution
      ),
      row.names = FALSE
    )
  },
  file = summary_file
)

session_information_file <- file.path(
  output_dir,
  "marker_supported_preservation_session_info.txt"
)

capture.output(
  sessionInfo(),
  file = session_information_file
)

output_files <- list.files(
  output_dir,
  full.names = TRUE
)

output_manifest <- data.frame(
  file_name = basename(output_files),
  file_size_bytes = file.info(
    output_files
  )$size,
  stringsAsFactors = FALSE
)

write.csv(
  output_manifest,
  file.path(
    output_dir,
    "marker_supported_preservation_output_manifest.csv"
  ),
  row.names = FALSE
)

required_final_files <- c(
  "Table_marker_supported_preservation_all_iterations.csv",
  "Table_marker_supported_preservation_summary.csv",
  "Table_marker_supported_population_metrics_all_iterations.csv",
  "Table_marker_supported_population_preservation_summary.csv",
  "marker_supported_cluster_mapping_all_iterations.csv",
  "cell_level_marker_supported_preservation_all_iterations.csv",
  "full_data_marker_supported_population_distribution.csv",
  "Figure_marker_supported_preservation_plot_data.csv",
  "Figure_marker_supported_preservation.png",
  "Figure_marker_supported_preservation.pdf",
  "marker_supported_preservation_input_checksums.csv",
  "marker_supported_preservation_summary.txt",
  "marker_supported_preservation_session_info.txt",
  "marker_supported_preservation_output_manifest.csv"
)

final_output_status <- file.exists(
  file.path(
    output_dir,
    required_final_files
  )
)

cat(
  "Required final outputs verified:",
  sum(final_output_status),
  "of",
  length(final_output_status),
  "\n"
)

stopifnot(all(final_output_status))