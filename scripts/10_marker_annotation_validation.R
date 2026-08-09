# ============================================================
# 10_marker_annotation_validation.R
#
# Purpose:
# Validate the resolution-0.6 cluster annotations using
# cluster marker genes and the existing reference labels.
# ============================================================

rm(list = ls())
gc()

options(stringsAsFactors = FALSE)
input_files <- c(
  labelled_object = "sce_labeled_majorityvote.rds",
  reference_object = "sce_final_clusters_locked_rpca_pcs20_res06.rds",
  marker_table = "markers_all_clusters_res06.csv"
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    "Missing required files: ",
    paste(missing_files, collapse = ", ")
  )
}

print(input_files)
cat("Loading labelled object...\n")
labelled_obj <- readRDS(
  unname(input_files["labelled_object"])
)

cat("Loading reference object...\n")
reference_obj <- readRDS(
  unname(input_files["reference_object"])
)

cat("Loading marker table...\n")
markers_res06 <- read.csv(
  unname(input_files["marker_table"]),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("Labelled object dimensions:", dim(labelled_obj), "\n")
cat("Reference object dimensions:", dim(reference_obj), "\n")
cat("Marker table dimensions:", dim(markers_res06), "\n")

output_dir <- file.path(
  "reviewer_revision_outputs",
  "10_marker_annotation_validation"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

cat("Output folder created:", output_dir, "\n")

stopifnot(
  identical(
    colnames(labelled_obj),
    colnames(reference_obj)
  )
)

stopifnot(
  "cluster_final" %in% colnames(labelled_obj@meta.data)
)

stopifnot(
  all(
    as.character(labelled_obj$cluster_final) ==
      as.character(reference_obj$cluster_final)
  )
)

cat(
  "Object alignment verified:",
  ncol(labelled_obj),
  "cells with identical cluster assignments.\n"
)
positive_markers <- markers_res06 %>%
  mutate(
    cluster = as.integer(cluster)
  ) %>%
  filter(
    !is.na(gene),
    p_val_adj < 0.05,
    avg_log2FC > 0.25,
    pct.1 > pct.2
  ) %>%
  arrange(
    cluster,
    desc(avg_log2FC)
  )

positive_marker_counts <- positive_markers %>%
  count(
    cluster,
    name = "number_of_positive_markers"
  )

print(positive_marker_counts)

top_marker_summary <- positive_markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 15,
    with_ties = FALSE
  ) %>%
  summarise(
    top_markers = paste(gene, collapse = ", "),
    .groups = "drop"
  ) %>%
  arrange(cluster)

print(
  as.data.frame(top_marker_summary),
  row.names = FALSE
)


marker_panels <- list(
  Endothelial = c(
    "PECAM1", "VWF", "KDR", "EMCN", "ENG",
    "CLDN5", "RAMP2", "PLVAP", "ESAM", "CDH5", "RGCC"
  ),
  Fibroblast_Stromal = c(
    "COL1A1", "COL1A2", "COL3A1", "DCN", "LUM",
    "PDGFRA", "THY1", "PRG4", "MMP1", "MMP3"
  ),
  Mural_Perivascular = c(
    "RGS5", "CSPG4", "MCAM", "PDGFRB", "RERGL",
    "ACTG2", "MYH11", "CNN1", "ABCC9", "GUCY1A2"
  ),
  Chondrocyte = c(
    "COL2A1", "ACAN", "SOX9", "COMP", "CHAD",
    "MATN3", "COL9A1", "COL9A2", "CRTAC1"
  ),
  T_Cell = c(
    "CD3D", "CD3E", "CD3G", "TRAC", "TRBC1",
    "TRBC2", "CD8A", "CD8B", "IL7R", "LTB", "NKG7"
  ),
  Myeloid = c(
    "LST1", "TYROBP", "FCER1G", "CTSS", "CD68",
    "C1QA", "C1QB", "C1QC", "AIF1", "APOE"
  ),
  Adipocyte = c(
    "ADIPOQ", "PLIN1", "FABP4", "LPL", "CIDEA",
    "CIDEC", "LEP", "PPARG", "DGAT2"
  ),
  Cycling = c(
    "MKI67", "TOP2A", "UBE2C", "AURKB",
    "CENPA", "NUSAP1", "HMMR"
  )
)

panel_availability <- data.frame(
  marker_panel = names(marker_panels),
  genes_requested = vapply(
    marker_panels,
    length,
    integer(1)
  ),
  genes_found_in_marker_table = vapply(
    marker_panels,
    function(x) {
      sum(x %in% markers_res06$gene)
    },
    integer(1)
  )
)

print(panel_availability, row.names = FALSE)


panel_gene_table <- bind_rows(
  lapply(
    names(marker_panels),
    function(panel_name) {
      data.frame(
        marker_panel = panel_name,
        gene = marker_panels[[panel_name]]
      )
    }
  )
)

panel_evidence <- positive_markers %>%
  inner_join(
    panel_gene_table,
    by = "gene"
  ) %>%
  group_by(
    cluster,
    marker_panel
  ) %>%
  summarise(
    positive_panel_genes = n_distinct(gene),
    evidence_genes = paste(
      sort(unique(gene)),
      collapse = ", "
    ),
    evidence_score = sum(
      avg_log2FC * (pct.1 - pct.2)
    ),
    .groups = "drop"
  )

all_cluster_panel_pairs <- expand_grid(
  cluster = sort(unique(markers_res06$cluster)),
  marker_panel = names(marker_panels)
)

panel_evidence_complete <- all_cluster_panel_pairs %>%
  left_join(
    panel_evidence,
    by = c("cluster", "marker_panel")
  ) %>%
  mutate(
    positive_panel_genes = replace_na(
      positive_panel_genes,
      0L
    ),
    evidence_genes = replace_na(
      evidence_genes,
      "None"
    ),
    evidence_score = replace_na(
      evidence_score,
      0
    )
  ) %>%
  arrange(
    cluster,
    desc(evidence_score)
  )

top_three_panels_per_cluster <- panel_evidence_complete %>%
  group_by(cluster) %>%
  slice_head(n = 3) %>%
  ungroup()

print(
  as.data.frame(top_three_panels_per_cluster),
  row.names = FALSE
)


genes_for_dotplot <- unique(
  unlist(
    marker_panels,
    use.names = FALSE
  )
)

genes_present <- genes_for_dotplot[
  genes_for_dotplot %in% rownames(labelled_obj)
]

genes_missing <- genes_for_dotplot[
  !genes_for_dotplot %in% rownames(labelled_obj)
]

cat(
  "Canonical genes present in object:",
  length(genes_present),
  "of",
  length(genes_for_dotplot),
  "\n"
)

cat(
  "Genes absent from object:",
  paste(genes_missing, collapse = ", "),
  "\n"
)


dotplot_panels <- list(
  Endothelial = c(
    "PECAM1", "VWF", "KDR", "EMCN", "CLDN5"
  ),
  Fibroblast_Stromal = c(
    "COL1A1", "COL1A2", "DCN", "LUM", "PRG4"
  ),
  Mural_Perivascular = c(
    "RGS5", "PDGFRB", "ACTG2", "MYH11", "ABCC9"
  ),
  Chondrocyte = c(
    "COL2A1", "ACAN", "SOX9", "COMP", "CRTAC1"
  ),
  T_Cell = c(
    "CD3D", "CD3E", "TRAC", "CD8A", "NKG7"
  ),
  Myeloid = c(
    "LST1", "TYROBP", "C1QA", "C1QC", "CD68"
  ),
  Adipocyte = c(
    "ADIPOQ", "PLIN1", "FABP4", "LPL", "PPARG"
  ),
  Cycling = c(
    "MKI67", "TOP2A", "UBE2C", "AURKB", "NUSAP1"
  )
)

dotplot_genes <- unique(
  unlist(
    dotplot_panels,
    use.names = FALSE
  )
)

stopifnot(
  all(dotplot_genes %in% rownames(labelled_obj))
)

cat(
  "DotPlot marker set verified:",
  length(dotplot_genes),
  "genes.\n"
)


labelled_obj$cluster_final <- factor(
  labelled_obj$cluster_final,
  levels = 0:17
)

Idents(labelled_obj) <- "cluster_final"

marker_dotplot <- DotPlot(
  object = labelled_obj,
  features = dotplot_panels,
  assay = "RNA",
  cols = c("grey90", "#B2182B"),
  dot.scale = 5
) +
  RotatedAxis() +
  labs(
    title = "Canonical marker expression across resolution-0.6 clusters",
    x = "Marker gene",
    y = "Cluster"
  ) +
  theme_bw(base_size = 10) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      size = 8
    ),
    axis.text.y = element_text(
      size = 9
    )
  )

print(marker_dotplot)

ggsave(
  filename = file.path(
    output_dir,
    "FigureS1_canonical_marker_dotplot.png"
  ),
  plot = marker_dotplot,
  width = 18,
  height = 8,
  dpi = 600,
  bg = "white",
  limitsize = FALSE
)

ggsave(
  filename = file.path(
    output_dir,
    "FigureS1_canonical_marker_dotplot.pdf"
  ),
  plot = marker_dotplot,
  width = 18,
  height = 8,
  bg = "white",
  limitsize = FALSE
)

cat("Marker DotPlot saved successfully.\n")


current_annotations <- labelled_obj@meta.data %>%
  transmute(
    cluster = as.integer(
      as.character(cluster_final)
    ),
    blueprint_label = label_blueprint_cluster,
    hpca_label = label_hpca_cluster,
    current_final_label = label_final
  ) %>%
  distinct() %>%
  arrange(cluster)

cluster_sizes <- labelled_obj@meta.data %>%
  transmute(
    cluster = as.integer(
      as.character(cluster_final)
    )
  ) %>%
  count(
    cluster,
    name = "number_of_cells"
  ) %>%
  arrange(cluster)

current_annotation_summary <- current_annotations %>%
  left_join(
    cluster_sizes,
    by = "cluster"
  )

print(
  as.data.frame(current_annotation_summary),
  row.names = FALSE
)

ranked_panel_summary <- top_three_panels_per_cluster %>%
  group_by(cluster) %>%
  mutate(
    evidence_rank = row_number()
  ) %>%
  ungroup() %>%
  select(
    cluster,
    evidence_rank,
    marker_panel,
    positive_panel_genes,
    evidence_genes,
    evidence_score
  ) %>%
  pivot_wider(
    names_from = evidence_rank,
    values_from = c(
      marker_panel,
      positive_panel_genes,
      evidence_genes,
      evidence_score
    ),
    names_glue = "{.value}_rank_{evidence_rank}"
  )

annotation_evidence_table <- current_annotation_summary %>%
  left_join(
    top_marker_summary,
    by = "cluster"
  ) %>%
  left_join(
    ranked_panel_summary,
    by = "cluster"
  ) %>%
  arrange(cluster)

write.csv(
  annotation_evidence_table,
  file.path(
    output_dir,
    "TableS1_preliminary_annotation_evidence.csv"
  ),
  row.names = FALSE
)

print(
  as.data.frame(
    annotation_evidence_table[
      ,
      c(
        "cluster",
        "current_final_label",
        "marker_panel_rank_1",
        "evidence_genes_rank_1"
      )
    ]
  ),
  row.names = FALSE
)



marker_supported_annotation <- data.frame(
  cluster = 0:17,
  marker_supported_label = c(
    "Fibroblast/stromal cells",  # 0
    "Endothelial cells",         # 1
    "Endothelial cells",         # 2
    "Endothelial cells",         # 3
    "Fibroblast/stromal cells",  # 4
    "Endothelial cells",         # 5
    "Endothelial cells",         # 6
    "Mural/perivascular cells",  # 7
    "Fibroblast/stromal cells",  # 8
    "T cells",                    # 9
    "Endothelial cells",         # 10
    "Myeloid cells",              # 11
    "Mural/perivascular cells",  # 12
    "Mural/perivascular cells",  # 13
    "Fibroblast/stromal cells",  # 14
    "T cells",                    # 15
    "Myeloid cells",              # 16
    "Myeloid cells"               # 17
  ),
  cluster_interpretation = c(
    "Matrix-remodelling fibroblast-like stromal cells",
    "Endothelial cells",
    "ACKR1-positive venular endothelial cells",
    "Activated endothelial cells",
    "PRG4/CRTAC1-positive lining stromal cells",
    "KDR/PLVAP-positive endothelial cells",
    "Inflammatory/activated endothelial cells",
    "Pericyte-like mural cells",
    "Matrix fibroblast-like stromal cells",
    "T-cell population with mixed T/NK-associated markers",
    "Arterial-like endothelial cells",
    "Monocyte/macrophage-like myeloid cells",
    "Smooth-muscle-like mural cells",
    "Perivascular/mural cells with endothelial overlap",
    "Cycling stromal cells",
    "T cells",
    "Macrophage-like myeloid cells",
    "Mixed myeloid/stromal cluster"
  ),
  annotation_confidence = c(
    "High", "High", "High", "High", "Moderate", "High",
    "High", "High", "High", "High", "High", "High",
    "High", "Moderate", "Moderate", "High", "High", "Low"
  ),
  marker_rationale = c(
    "COL1A1, COL1A2, COL3A1, MMP1, MMP3 and PRG4",
    "PECAM1, VWF, EMCN, ENG, ESAM and CDH5",
    "ACKR1, VWF, PECAM1 and PLVAP",
    "PECAM1, VWF, EMCN, PLVAP, KDR and ENG",
    "PRG4, CRTAC1, COL1A1, COL1A2, COL3A1 and PDGFRA",
    "KDR, PLVAP, CLDN5, EMCN, CDH5 and RAMP2",
    "PECAM1, VWF, PLVAP, EMCN, SELE and ENG",
    "RGS5, PDGFRB, MCAM, ABCC9, GUCY1A2 and MYH11",
    "DCN, LUM, COL1A1, COL1A2, COL3A1 and PDGFRA",
    "CD3D, CD3E, TRAC, TRBC1, TRBC2, IL7R and NKG7",
    "GJA5, CLDN5, PECAM1, KDR, CDH5 and RAMP2",
    "LST1, TYROBP, C1QA, C1QB, C1QC, CD68 and FCER1G",
    "MYH11, ACTG2, CNN1, RERGL, PDGFRB and RGS5",
    "RGS5, ABCC9, GUCY1A2, PDGFRB with endothelial markers",
    "MKI67, TOP2A, UBE2C, AURKB and stromal collagens",
    "CD3D, CD3E, CD3G, TRBC2, IL7R and LTB",
    "C1QA, C1QB, C1QC, TYROBP, LST1, CD68 and APOE",
    "Myeloid markers with concurrent PRG4/DCN/LUM expression"
  ),
  stringsAsFactors = FALSE
)

annotation_review_table <- annotation_evidence_table %>%
  left_join(
    marker_supported_annotation,
    by = "cluster"
  ) %>%
  arrange(cluster)

print(
  annotation_review_table[
    ,
    c(
      "cluster",
      "current_final_label",
      "marker_supported_label",
      "annotation_confidence"
    )
  ],
  row.names = FALSE
)

annotation_review_file <- file.path(
  output_dir,
  "TableS2_marker_supported_annotation_review.csv"
)

write.csv(
  annotation_review_table,
  annotation_review_file,
  row.names = FALSE
)

cat(
  "Annotation-review table saved:",
  annotation_review_file,
  "\nFile exists:",
  file.exists(annotation_review_file),
  "\n"
)

labelled_obj_marker_validated <- labelled_obj

annotation_map <- setNames(
  marker_supported_annotation$marker_supported_label,
  as.character(marker_supported_annotation$cluster)
)

interpretation_map <- setNames(
  marker_supported_annotation$cluster_interpretation,
  as.character(marker_supported_annotation$cluster)
)

confidence_map <- setNames(
  marker_supported_annotation$annotation_confidence,
  as.character(marker_supported_annotation$cluster)
)

cell_clusters <- as.character(
  labelled_obj_marker_validated$cluster_final
)

stopifnot(
  all(cell_clusters %in% names(annotation_map))
)

labelled_obj_marker_validated$marker_supported_label <-
  unname(annotation_map[cell_clusters])

labelled_obj_marker_validated$marker_interpretation <-
  unname(interpretation_map[cell_clusters])

labelled_obj_marker_validated$annotation_confidence <-
  unname(confidence_map[cell_clusters])

stopifnot(
  !anyNA(labelled_obj_marker_validated$marker_supported_label),
  !anyNA(labelled_obj_marker_validated$marker_interpretation),
  !anyNA(labelled_obj_marker_validated$annotation_confidence)
)

print(
  table(
    labelled_obj_marker_validated$marker_supported_label,
    useNA = "ifany"
  )
)

cat("Marker-supported annotations added to the copied object.\n")


validated_object_file <- file.path(
  output_dir,
  "sce_marker_validated_broad_annotations_res06.rds"
)

saveRDS(
  labelled_obj_marker_validated,
  validated_object_file,
  compress = "gzip"
)

cat(
  "Validated object saved:",
  validated_object_file,
  "\nFile exists:",
  file.exists(validated_object_file),
  "\nFile size (GB):",
  round(file.info(validated_object_file)$size / 1024^3, 3),
  "\n"
)



if (!requireNamespace("patchwork", quietly = TRUE)) {
  stop("Install patchwork using: install.packages('patchwork')")
}

library(patchwork)

stopifnot(
  "umap.rpca" %in% Reductions(labelled_obj_marker_validated)
)

labelled_obj_marker_validated$marker_supported_label <- factor(
  labelled_obj_marker_validated$marker_supported_label,
  levels = c(
    "Endothelial cells",
    "Fibroblast/stromal cells",
    "Mural/perivascular cells",
    "T cells",
    "Myeloid cells"
  )
)

annotation_colours <- c(
  "Endothelial cells" = "#0072B2",
  "Fibroblast/stromal cells" = "#E69F00",
  "Mural/perivascular cells" = "#009E73",
  "T cells" = "#CC79A7",
  "Myeloid cells" = "#D55E00"
)

umap_clusters <- DimPlot(
  labelled_obj_marker_validated,
  reduction = "umap.rpca",
  group.by = "cluster_final",
  label = TRUE,
  repel = TRUE,
  label.size = 4,
  pt.size = 0.15,
  raster = FALSE
) +
  NoLegend() +
  labs(
    title = "A. Resolution-0.6 clusters",
    x = "UMAP 1",
    y = "UMAP 2"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13)
  )

umap_annotations <- DimPlot(
  labelled_obj_marker_validated,
  reduction = "umap.rpca",
  group.by = "marker_supported_label",
  cols = annotation_colours,
  pt.size = 0.15,
  raster = FALSE
) +
  labs(
    title = "B. Marker-supported broad annotations",
    x = "UMAP 1",
    y = "UMAP 2",
    colour = "Cell population"
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    legend.title = element_text(face = "bold")
  )

marker_validated_umap <- umap_clusters + umap_annotations

print(marker_validated_umap)

ggsave(
  file.path(output_dir, "FigureS2_marker_validated_UMAP.png"),
  marker_validated_umap,
  width = 14,
  height = 6.5,
  dpi = 400,
  bg = "white"
)

ggsave(
  file.path(output_dir, "FigureS2_marker_validated_UMAP.pdf"),
  marker_validated_umap,
  width = 14,
  height = 6.5,
  bg = "white"
)

cat("Marker-supported UMAP saved successfully.\n")



browseURL(
  normalizePath(
    file.path(
      output_dir,
      "FigureS2_marker_validated_UMAP.png"
    )
  )
)


annotation_population_summary <-
  labelled_obj_marker_validated@meta.data %>%
  count(
    marker_supported_label,
    name = "number_of_cells",
    .drop = FALSE
  ) %>%
  mutate(
    percentage_of_total_cells = round(
      100 * number_of_cells / sum(number_of_cells),
      2
    )
  ) %>%
  arrange(desc(number_of_cells))

stopifnot(
  sum(annotation_population_summary$number_of_cells) == 43110
)

write.csv(
  annotation_population_summary,
  file.path(
    output_dir,
    "TableS3_marker_supported_population_summary.csv"
  ),
  row.names = FALSE
)

print(
  as.data.frame(annotation_population_summary),
  row.names = FALSE
)

cat("Population summary saved successfully.\n")


donor_annotation_summary <-
  labelled_obj_marker_validated@meta.data %>%
  count(
    sample_id,
    marker_supported_label,
    name = "number_of_cells",
    .drop = FALSE
  ) %>%
  group_by(sample_id) %>%
  mutate(
    total_cells_in_donor = sum(number_of_cells),
    percentage_within_donor = round(
      100 * number_of_cells / total_cells_in_donor,
      2
    )
  ) %>%
  ungroup() %>%
  arrange(sample_id, marker_supported_label)

stopifnot(
  sum(donor_annotation_summary$number_of_cells) == 43110
)

write.csv(
  donor_annotation_summary,
  file.path(
    output_dir,
    "TableS4_marker_supported_populations_by_donor.csv"
  ),
  row.names = FALSE
)

print(
  as.data.frame(donor_annotation_summary),
  row.names = FALSE
)

cat("Donor-level annotation summary saved successfully.\n")



annotation_confidence_summary <-
  labelled_obj_marker_validated@meta.data %>%
  group_by(annotation_confidence) %>%
  summarise(
    number_of_clusters = n_distinct(cluster_final),
    number_of_cells = n(),
    percentage_of_total_cells = round(
      100 * n() / nrow(labelled_obj_marker_validated@meta.data),
      2
    ),
    .groups = "drop"
  ) %>%
  mutate(
    annotation_confidence = factor(
      annotation_confidence,
      levels = c("High", "Moderate", "Low")
    )
  ) %>%
  arrange(annotation_confidence)

write.csv(
  annotation_confidence_summary,
  file.path(
    output_dir,
    "TableS5_annotation_confidence_summary.csv"
  ),
  row.names = FALSE
)

print(
  as.data.frame(annotation_confidence_summary),
  row.names = FALSE
)

cat("Annotation-confidence summary saved successfully.\n")

session_info_file <- file.path(
  output_dir,
  "marker_annotation_validation_session_info.txt"
)

capture.output(
  sessionInfo(),
  file = session_info_file
)

cat(
  "Session information saved:",
  session_info_file,
  "\nFile exists:",
  file.exists(session_info_file),
  "\n"
)



