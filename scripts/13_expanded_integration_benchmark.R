# ============================================================
# 13_expanded_integration_benchmark.R
#
# Purpose:
# Compare unintegrated PCA, CCA, RPCA, and Harmony using the
# same donor-balanced cells and common preprocessing settings.
# Metrics will include clustering agreement, silhouette,
# kBET, iLISI, graph connectivity, and runtime.
# ============================================================


# ============================================================
# 1. Clean session and load required packages
# ============================================================

rm(list = ls())
gc()

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "mclust",
  "aricode",
  "cluster",
  "dplyr",
  "tidyr",
  "ggplot2",
  "future",
  "harmony",
  "RANN",
  "igraph"
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
  library(cluster)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(future)
  library(harmony)
})

future::plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

benchmark_seed <- 2026
set.seed(benchmark_seed)

package_versions <- data.frame(
  package = required_packages,
  version = vapply(
    required_packages,
    function(package_name) {
      as.character(packageVersion(package_name))
    },
    character(1)
  ),
  stringsAsFactors = FALSE
)

print(package_versions, row.names = FALSE)
cat("Expanded-benchmark packages loaded successfully.\n")





# ============================================================
# 2. Define and verify input/output files
# ============================================================

output_dir <- file.path(
  "reviewer_revision_outputs",
  "13_expanded_integration_benchmark"
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
    "sce_final_clusters_locked_rpca_pcs20_res06.rds",
  marker_validated_object = file.path(
    "reviewer_revision_outputs",
    "10_marker_annotation_validation",
    "sce_marker_validated_broad_annotations_res06.rds"
  ),
  benchmark_cells = file.path(
    "workflow_benchmark_outputs",
    "benchmark_sampled_cells.csv"
  )
)

input_check <- data.frame(
  input_name = names(input_files),
  file_path = unname(input_files),
  exists = file.exists(input_files),
  stringsAsFactors = FALSE
)

print(input_check, row.names = FALSE)

if (any(!input_check$exists)) {
  stop(
    "Missing required files: ",
    paste(
      input_check$input_name[!input_check$exists],
      collapse = ", "
    )
  )
}

cat("All expanded-benchmark input files verified.\n")
cat("Output folder:", output_dir, "\n")



# ============================================================
# 3. Load objects and verify the common benchmark cells
# ============================================================

cat("Loading post-filter object...\n")

post_obj <- readRDS(
  unname(input_files["post_filter_object"])
)

stopifnot(
  inherits(post_obj, "Seurat"),
  "sample_id" %in% colnames(post_obj@meta.data)
)

DefaultAssay(post_obj) <- "RNA"

cat("Loading reference clustering...\n")

reference_temp <- readRDS(
  unname(input_files["reference_object"])
)

stopifnot(
  inherits(reference_temp, "Seurat"),
  identical(colnames(post_obj), colnames(reference_temp)),
  "cluster_final" %in% colnames(reference_temp@meta.data)
)

reference_labels <- as.character(
  reference_temp@meta.data$cluster_final
)

names(reference_labels) <- rownames(
  reference_temp@meta.data
)

rm(reference_temp)
invisible(gc())

cat("Loading marker-supported annotations...\n")

validated_temp <- readRDS(
  unname(input_files["marker_validated_object"])
)

stopifnot(
  inherits(validated_temp, "Seurat"),
  identical(colnames(post_obj), colnames(validated_temp)),
  "marker_supported_label" %in%
    colnames(validated_temp@meta.data)
)

biological_labels <- as.character(
  validated_temp@meta.data$marker_supported_label
)

names(biological_labels) <- rownames(
  validated_temp@meta.data
)

rm(validated_temp)
invisible(gc())

benchmark_map <- read.csv(
  unname(input_files["benchmark_cells"]),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

required_mapping_columns <- c(
  "internal_cell_id",
  "original_cell_barcode",
  "donor_id"
)

missing_mapping_columns <- setdiff(
  required_mapping_columns,
  colnames(benchmark_map)
)

if (length(missing_mapping_columns) > 0) {
  stop(
    "Missing benchmark-map columns: ",
    paste(missing_mapping_columns, collapse = ", ")
  )
}

if (
  anyDuplicated(benchmark_map$original_cell_barcode) > 0 ||
  anyDuplicated(benchmark_map$internal_cell_id) > 0
) {
  stop("Duplicate cells detected in the benchmark map.")
}

if (
  !all(
    benchmark_map$original_cell_barcode %in%
    colnames(post_obj)
  )
) {
  stop("Some benchmark cells are absent from post_obj.")
}

post_donors <- as.character(
  post_obj@meta.data$sample_id
)

names(post_donors) <- rownames(
  post_obj@meta.data
)

mapped_donors <- post_donors[
  benchmark_map$original_cell_barcode
]

if (
  !identical(
    unname(mapped_donors),
    benchmark_map$donor_id
  )
) {
  stop("Donor labels do not match the benchmark map.")
}

benchmark_donor_counts <- table(
  benchmark_map$donor_id
)

stopifnot(
  nrow(benchmark_map) == 6955,
  length(benchmark_donor_counts) == 5,
  length(unique(as.integer(benchmark_donor_counts))) == 1
)

cat(
  "Source object dimensions:",
  nrow(post_obj), "genes and",
  ncol(post_obj), "cells.\n"
)

cat(
  "Benchmark subset verified:",
  nrow(benchmark_map),
  "cells.\n"
)

print(benchmark_donor_counts)

cat("Reference and biological labels verified.\n")





# ============================================================
# 4. Construct the aligned benchmark count matrix
# ============================================================

benchmark_cells <- benchmark_map$original_cell_barcode
benchmark_internal_ids <- benchmark_map$internal_cell_id
feature_universe <- rownames(post_obj)

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
    "Expected five count layers, but found ",
    length(count_layer_names), "."
  )
}

extract_benchmark_count_layer <- function(layer_name) {
  
  cat("Extracting:", layer_name, "\n")
  
  layer_matrix <- LayerData(
    object = post_obj,
    assay = "RNA",
    layer = layer_name
  )
  
  selected_layer_cells <- colnames(layer_matrix)[
    colnames(layer_matrix) %in% benchmark_cells
  ]
  
  if (length(selected_layer_cells) == 0) {
    return(NULL)
  }
  
  unexpected_features <- setdiff(
    rownames(layer_matrix),
    feature_universe
  )
  
  if (length(unexpected_features) > 0) {
    stop(
      "Unexpected features detected in ",
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
      ncol = length(selected_layer_cells),
      sparse = TRUE,
      dimnames = list(
        missing_features,
        selected_layer_cells
      )
    )
    
    layer_matrix <- rbind(
      layer_matrix[, selected_layer_cells, drop = FALSE],
      zero_matrix
    )
    
  } else {
    
    layer_matrix <- layer_matrix[
      ,
      selected_layer_cells,
      drop = FALSE
    ]
  }
  
  layer_matrix[
    feature_universe,
    ,
    drop = FALSE
  ]
}

benchmark_count_list <- lapply(
  count_layer_names,
  extract_benchmark_count_layer
)

benchmark_count_list <- Filter(
  Negate(is.null),
  benchmark_count_list
)

benchmark_counts <- do.call(
  cbind,
  benchmark_count_list
)

if (anyDuplicated(colnames(benchmark_counts)) > 0) {
  stop("Duplicate cell barcodes detected.")
}

if (!setequal(colnames(benchmark_counts), benchmark_cells)) {
  stop("The reconstructed cell set does not match the benchmark.")
}

benchmark_counts <- benchmark_counts[
  ,
  benchmark_cells,
  drop = FALSE
]

colnames(benchmark_counts) <- benchmark_internal_ids

benchmark_metadata <- post_obj@meta.data[
  benchmark_cells,
  ,
  drop = FALSE
]

benchmark_metadata$original_cell_barcode <-
  benchmark_cells

benchmark_metadata$donor_id_analysis <-
  benchmark_map$donor_id

benchmark_metadata$reference_cluster <-
  unname(reference_labels[benchmark_cells])

benchmark_metadata$biological_label <-
  unname(biological_labels[benchmark_cells])

rownames(benchmark_metadata) <-
  benchmark_internal_ids

stopifnot(
  identical(
    colnames(benchmark_counts),
    rownames(benchmark_metadata)
  ),
  !anyNA(benchmark_metadata$reference_cluster),
  !anyNA(benchmark_metadata$biological_label)
)

write.csv(
  data.frame(
    internal_cell_id = benchmark_internal_ids,
    original_cell_barcode = benchmark_cells,
    donor_id = benchmark_metadata$donor_id_analysis,
    reference_cluster = benchmark_metadata$reference_cluster,
    biological_label = benchmark_metadata$biological_label
  ),
  file.path(
    output_dir,
    "benchmark_cells_with_reference_labels.csv"
  ),
  row.names = FALSE
)

rm(
  post_obj,
  benchmark_count_list
)

invisible(gc())

cat(
  "Benchmark count matrix:",
  paste(dim(benchmark_counts), collapse = " x "),
  "\n"
)

cat("Reference clusters represented:\n")
print(table(benchmark_metadata$reference_cluster))

cat("Marker-supported populations represented:\n")
print(table(benchmark_metadata$biological_label))

cat("Benchmark matrix and metadata verified.\n")

# ============================================================
# 5. Define common workflow functions
# ============================================================

create_preprocessed_object <- function(
    split_by_donor = FALSE
) {
  
  obj <- CreateSeuratObject(
    counts = benchmark_counts,
    meta.data = benchmark_metadata,
    assay = "RNA",
    project = "ExpandedIntegrationBenchmark",
    min.cells = 0,
    min.features = 0
  )
  
  DefaultAssay(obj) <- "RNA"
  
  if (split_by_donor) {
    
    obj[["RNA"]] <- split(
      obj[["RNA"]],
      f = obj$donor_id_analysis
    )
  }
  
  set.seed(benchmark_seed)
  
  obj <- NormalizeData(
    obj,
    normalization.method = "LogNormalize",
    scale.factor = 10000,
    verbose = FALSE
  )
  
  obj <- FindVariableFeatures(
    obj,
    selection.method = "vst",
    nfeatures = 2000,
    verbose = FALSE
  )
  
  obj <- ScaleData(
    obj,
    features = VariableFeatures(obj),
    verbose = FALSE
  )
  
  obj <- RunPCA(
    obj,
    features = VariableFeatures(obj),
    npcs = 20,
    seed.use = benchmark_seed,
    verbose = FALSE
  )
  
  obj
}


cluster_from_reduction <- function(
    obj,
    reduction_name
) {
  
  set.seed(benchmark_seed)
  
  obj <- FindNeighbors(
    obj,
    reduction = reduction_name,
    dims = 1:20,
    k.param = 20,
    verbose = FALSE
  )
  
  obj <- FindClusters(
    obj,
    resolution = 0.6,
    algorithm = 3,
    random.seed = benchmark_seed,
    verbose = FALSE
  )
  
  cluster_labels <- as.character(
    Idents(obj)
  )
  
  names(cluster_labels) <- colnames(obj)
  
  embedding_matrix <- Embeddings(
    obj,
    reduction = reduction_name
  )[
    ,
    1:20,
    drop = FALSE
  ]
  
  list(
    labels = cluster_labels,
    embedding = embedding_matrix,
    clusters_recovered =
      length(unique(cluster_labels))
  )
}


safe_nmi <- function(
    labels_a,
    labels_b
) {
  
  labels_a <- as.character(labels_a)
  labels_b <- as.character(labels_b)
  
  combined_levels <- sort(
    unique(
      c(labels_a, labels_b)
    )
  )
  
  labels_a_integer <- as.integer(
    factor(
      labels_a,
      levels = combined_levels
    )
  )
  
  labels_b_integer <- as.integer(
    factor(
      labels_b,
      levels = combined_levels
    )
  )
  
  aricode::NMI(
    labels_a_integer,
    labels_b_integer
  )
}

cat("Common preprocessing and clustering functions defined.\n")


# ============================================================
# 6A. Run the unintegrated PCA workflow
# ============================================================

run_unintegrated_workflow <- function() {
  
  start_time <- Sys.time()
  
  obj <- create_preprocessed_object(
    split_by_donor = FALSE
  )
  
  clustered_result <- cluster_from_reduction(
    obj = obj,
    reduction_name = "pca"
  )
  
  runtime_seconds <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )
  
  list(
    workflow = "Unintegrated PCA",
    labels = clustered_result$labels,
    embedding = clustered_result$embedding,
    clusters_recovered =
      clustered_result$clusters_recovered,
    runtime_seconds = runtime_seconds
  )
}

cat("Running unintegrated PCA workflow...\n")

unintegrated_result <-
  run_unintegrated_workflow()

saveRDS(
  unintegrated_result,
  file.path(
    output_dir,
    "unintegrated_PCA_result.rds"
  )
)

print(
  data.frame(
    workflow = unintegrated_result$workflow,
    clusters_recovered =
      unintegrated_result$clusters_recovered,
    runtime_seconds = round(
      unintegrated_result$runtime_seconds,
      2
    )
  ),
  row.names = FALSE
)

cat("Unintegrated PCA workflow completed and saved.\n")



# ============================================================
# 6B. Run the CCA-integrated workflow
# ============================================================

run_cca_workflow <- function() {
  
  start_time <- Sys.time()
  
  obj <- create_preprocessed_object(
    split_by_donor = TRUE
  )
  
  set.seed(benchmark_seed)
  
  obj <- IntegrateLayers(
    object = obj,
    method = CCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.cca",
    dims = 1:20,
    verbose = FALSE
  )
  
  clustered_result <- cluster_from_reduction(
    obj = obj,
    reduction_name = "integrated.cca"
  )
  
  runtime_seconds <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )
  
  list(
    workflow = "CCA",
    labels = clustered_result$labels,
    embedding = clustered_result$embedding,
    clusters_recovered =
      clustered_result$clusters_recovered,
    runtime_seconds = runtime_seconds
  )
}

cat("Running CCA-integrated workflow...\n")

cca_result <- run_cca_workflow()

saveRDS(
  cca_result,
  file.path(
    output_dir,
    "CCA_integration_result.rds"
  )
)

print(
  data.frame(
    workflow = cca_result$workflow,
    clusters_recovered =
      cca_result$clusters_recovered,
    runtime_seconds = round(
      cca_result$runtime_seconds,
      2
    )
  ),
  row.names = FALSE
)

cat("CCA-integrated workflow completed and saved.\n")



# ============================================================
# 6C. Run the RPCA-integrated workflow
# ============================================================

run_rpca_workflow <- function() {
  
  start_time <- Sys.time()
  
  obj <- create_preprocessed_object(
    split_by_donor = TRUE
  )
  
  set.seed(benchmark_seed)
  
  obj <- IntegrateLayers(
    object = obj,
    method = RPCAIntegration,
    orig.reduction = "pca",
    new.reduction = "integrated.rpca",
    dims = 1:20,
    verbose = FALSE
  )
  
  clustered_result <- cluster_from_reduction(
    obj = obj,
    reduction_name = "integrated.rpca"
  )
  
  runtime_seconds <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )
  
  list(
    workflow = "RPCA",
    labels = clustered_result$labels,
    embedding = clustered_result$embedding,
    clusters_recovered =
      clustered_result$clusters_recovered,
    runtime_seconds = runtime_seconds
  )
}

cat("Running RPCA-integrated workflow...\n")

rpca_result <- run_rpca_workflow()

saveRDS(
  rpca_result,
  file.path(
    output_dir,
    "RPCA_integration_result.rds"
  )
)

print(
  data.frame(
    workflow = rpca_result$workflow,
    clusters_recovered =
      rpca_result$clusters_recovered,
    runtime_seconds = round(
      rpca_result$runtime_seconds,
      2
    )
  ),
  row.names = FALSE
)

cat("RPCA-integrated workflow completed and saved.\n")



# ============================================================
# 6D. Run the Harmony-integrated workflow
# ============================================================

run_harmony_workflow <- function() {
  
  start_time <- Sys.time()
  
  obj <- create_preprocessed_object(
    split_by_donor = FALSE
  )
  
  set.seed(benchmark_seed)
  
  obj <- harmony::RunHarmony(
    object = obj,
    group.by.vars = "donor_id_analysis",
    reduction.use = "pca",
    dims.use = 1:20,
    reduction.save = "harmony",
    project.dim = FALSE
  )
  
  if (!"harmony" %in% Reductions(obj)) {
    stop("Harmony reduction was not created.")
  }
  
  clustered_result <- cluster_from_reduction(
    obj = obj,
    reduction_name = "harmony"
  )
  
  runtime_seconds <- as.numeric(
    difftime(
      Sys.time(),
      start_time,
      units = "secs"
    )
  )
  
  list(
    workflow = "Harmony",
    labels = clustered_result$labels,
    embedding = clustered_result$embedding,
    clusters_recovered =
      clustered_result$clusters_recovered,
    runtime_seconds = runtime_seconds
  )
}

cat("Running Harmony-integrated workflow...\n")

harmony_result <- run_harmony_workflow()

saveRDS(
  harmony_result,
  file.path(
    output_dir,
    "Harmony_integration_result.rds"
  )
)

print(
  data.frame(
    workflow = harmony_result$workflow,
    clusters_recovered =
      harmony_result$clusters_recovered,
    runtime_seconds = round(
      harmony_result$runtime_seconds,
      2
    )
  ),
  row.names = FALSE
)

cat("Harmony-integrated workflow completed and saved.\n")




# ============================================================
# 7. Verify workflow alignment and calculate initial metrics
# ============================================================

workflow_results <- list(
  "Unintegrated PCA" = unintegrated_result,
  "CCA" = cca_result,
  "RPCA" = rpca_result,
  "Harmony" = harmony_result
)

expected_cell_ids <- rownames(
  benchmark_metadata
)

for (workflow_name in names(workflow_results)) {
  
  current_result <- getElement(
    workflow_results,
    workflow_name
  )
  
  if (
    !setequal(
      names(current_result$labels),
      expected_cell_ids
    )
  ) {
    stop(
      "Cluster-label mismatch for ",
      workflow_name
    )
  }
  
  if (
    !setequal(
      rownames(current_result$embedding),
      expected_cell_ids
    )
  ) {
    stop(
      "Embedding mismatch for ",
      workflow_name
    )
  }
  
  current_result$labels <-
    current_result$labels[
      expected_cell_ids
    ]
  
  current_result$embedding <-
    current_result$embedding[
      expected_cell_ids,
      ,
      drop = FALSE
    ]
  
  if (
    anyNA(current_result$labels) ||
    anyNA(current_result$embedding)
  ) {
    stop(
      "Missing workflow values for ",
      workflow_name
    )
  }
  
  workflow_results[workflow_name] <-
    list(current_result)
}

reference_cluster_vector <- as.character(
  benchmark_metadata$reference_cluster
)

biological_label_vector <- as.character(
  benchmark_metadata$biological_label
)

set.seed(benchmark_seed)

silhouette_cell_ids <- sample(
  expected_cell_ids,
  size = min(
    3000,
    length(expected_cell_ids)
  ),
  replace = FALSE
)

calculate_mean_silhouette <- function(
    embedding_matrix,
    label_vector,
    selected_cell_ids
) {
  
  selected_embedding <- embedding_matrix[
    selected_cell_ids,
    ,
    drop = FALSE
  ]
  
  selected_labels <- label_vector[
    selected_cell_ids
  ]
  
  label_integer <- as.integer(
    factor(selected_labels)
  )
  
  silhouette_result <- cluster::silhouette(
    label_integer,
    dist(selected_embedding)
  )
  
  c(
    mean = mean(
      silhouette_result[, "sil_width"]
    ),
    median = median(
      silhouette_result[, "sil_width"]
    )
  )
}

initial_metric_list <- lapply(
  names(workflow_results),
  function(workflow_name) {
    
    current_result <- getElement(
      workflow_results,
      workflow_name
    )
    
    workflow_labels <- current_result$labels
    
    cluster_silhouette <- calculate_mean_silhouette(
      embedding_matrix = current_result$embedding,
      label_vector = workflow_labels,
      selected_cell_ids = silhouette_cell_ids
    )
    
    biological_silhouette <- calculate_mean_silhouette(
      embedding_matrix = current_result$embedding,
      label_vector = setNames(
        biological_label_vector,
        expected_cell_ids
      ),
      selected_cell_ids = silhouette_cell_ids
    )
    
    data.frame(
      workflow = workflow_name,
      clusters_recovered =
        current_result$clusters_recovered,
      ARI_vs_full_RPCA_reference =
        mclust::adjustedRandIndex(
          workflow_labels,
          reference_cluster_vector
        ),
      NMI_vs_full_RPCA_reference =
        safe_nmi(
          workflow_labels,
          reference_cluster_vector
        ),
      mean_cluster_silhouette =
        unname(
          cluster_silhouette["mean"]
        ),
      mean_biological_silhouette =
        unname(
          biological_silhouette["mean"]
        ),
      runtime_seconds =
        current_result$runtime_seconds,
      stringsAsFactors = FALSE
    )
  }
)

initial_metrics <- bind_rows(
  initial_metric_list
)

write.csv(
  initial_metrics,
  file.path(
    output_dir,
    "Table_initial_workflow_metrics.csv"
  ),
  row.names = FALSE
)

cell_level_assignments <- data.frame(
  internal_cell_id = expected_cell_ids,
  original_cell_barcode =
    benchmark_metadata$original_cell_barcode,
  donor_id =
    benchmark_metadata$donor_id_analysis,
  reference_cluster =
    reference_cluster_vector,
  biological_label =
    biological_label_vector,
  unintegrated_PCA_cluster =
    getElement(
      workflow_results,
      "Unintegrated PCA"
    )$labels,
  CCA_cluster =
    getElement(
      workflow_results,
      "CCA"
    )$labels,
  RPCA_cluster =
    getElement(
      workflow_results,
      "RPCA"
    )$labels,
  Harmony_cluster =
    getElement(
      workflow_results,
      "Harmony"
    )$labels,
  stringsAsFactors = FALSE
)

write.csv(
  cell_level_assignments,
  file.path(
    output_dir,
    "cell_level_assignments_all_workflows.csv"
  ),
  row.names = FALSE
)

print(
  initial_metrics %>%
    mutate(
      across(
        where(is.numeric),
        function(x) round(x, 4)
      )
    ),
  row.names = FALSE
)

cat(
  "Initial agreement and silhouette metrics saved successfully.\n"
)





# ============================================================
# 8A. Define kNN, iLISI, kBET, and graph-connectivity functions
# ============================================================

build_knn_result <- function(
    embedding_matrix,
    requested_neighbors = 90
) {
  
  if (nrow(embedding_matrix) <= requested_neighbors) {
    stop("Insufficient cells for the requested neighbourhood.")
  }
  
  raw_knn <- RANN::nn2(
    data = embedding_matrix,
    query = embedding_matrix,
    k = requested_neighbors,
    eps = 0
  )
  
  retained_neighbors <- requested_neighbors - 1
  
  neighbor_indices <- matrix(
    NA_integer_,
    nrow = nrow(embedding_matrix),
    ncol = retained_neighbors
  )
  
  neighbor_distances <- matrix(
    NA_real_,
    nrow = nrow(embedding_matrix),
    ncol = retained_neighbors
  )
  
  for (cell_index in seq_len(nrow(embedding_matrix))) {
    
    non_self_positions <- which(
      raw_knn$nn.idx[cell_index, ] != cell_index
    )
    
    if (
      length(non_self_positions) <
      retained_neighbors
    ) {
      stop(
        "Unable to recover sufficient non-self neighbours."
      )
    }
    
    retained_positions <- non_self_positions[
      seq_len(retained_neighbors)
    ]
    
    neighbor_indices[cell_index, ] <-
      raw_knn$nn.idx[
        cell_index,
        retained_positions
      ]
    
    neighbor_distances[cell_index, ] <-
      raw_knn$nn.dists[
        cell_index,
        retained_positions
      ]
  }
  
  rownames(neighbor_indices) <-
    rownames(embedding_matrix)
  
  rownames(neighbor_distances) <-
    rownames(embedding_matrix)
  
  list(
    indices = neighbor_indices,
    distances = neighbor_distances
  )
}


calculate_lisi_weights <- function(
    distance_vector,
    perplexity = 30,
    tolerance = 1e-5,
    maximum_iterations = 50
) {
  
  target_entropy <- log(perplexity)
  
  beta <- 1
  beta_minimum <- -Inf
  beta_maximum <- Inf
  
  calculate_entropy <- function(current_beta) {
    
    unnormalised_probability <- exp(
      -distance_vector * current_beta
    )
    
    probability_sum <- sum(
      unnormalised_probability
    )
    
    if (
      probability_sum == 0 ||
      !is.finite(probability_sum)
    ) {
      return(
        list(
          entropy = 0,
          probabilities = rep(
            NA_real_,
            length(distance_vector)
          )
        )
      )
    }
    
    entropy <- log(probability_sum) +
      current_beta *
      sum(
        distance_vector *
          unnormalised_probability
      ) /
      probability_sum
    
    list(
      entropy = entropy,
      probabilities =
        unnormalised_probability /
        probability_sum
    )
  }
  
  entropy_result <- calculate_entropy(beta)
  
  for (
    iteration_number in
    seq_len(maximum_iterations)
  ) {
    
    entropy_difference <-
      entropy_result$entropy -
      target_entropy
    
    if (
      !is.finite(entropy_difference) ||
      abs(entropy_difference) <= tolerance
    ) {
      break
    }
    
    if (entropy_difference > 0) {
      
      beta_minimum <- beta
      
      if (!is.finite(beta_maximum)) {
        beta <- beta * 2
      } else {
        beta <- (
          beta + beta_maximum
        ) / 2
      }
      
    } else {
      
      beta_maximum <- beta
      
      if (!is.finite(beta_minimum)) {
        beta <- beta / 2
      } else {
        beta <- (
          beta + beta_minimum
        ) / 2
      }
    }
    
    entropy_result <- calculate_entropy(beta)
  }
  
  entropy_result$probabilities
}


calculate_ilisi <- function(
    knn_result,
    batch_labels,
    perplexity = 30
) {
  
  batch_factor <- factor(
    batch_labels
  )
  
  batch_codes <- as.integer(
    batch_factor
  )
  
  number_of_batches <- nlevels(
    batch_factor
  )
  
  ilisi_values <- numeric(
    nrow(knn_result$indices)
  )
  
  for (
    cell_index in
    seq_len(nrow(knn_result$indices))
  ) {
    
    local_probabilities <- calculate_lisi_weights(
      distance_vector =
        knn_result$distances[cell_index, ],
      perplexity = perplexity
    )
    
    if (anyNA(local_probabilities)) {
      ilisi_values[cell_index] <- NA_real_
      next
    }
    
    local_batch_codes <- batch_codes[
      knn_result$indices[cell_index, ]
    ]
    
    batch_probabilities <- vapply(
      seq_len(number_of_batches),
      function(batch_number) {
        sum(
          local_probabilities[
            local_batch_codes ==
              batch_number
          ]
        )
      },
      numeric(1)
    )
    
    ilisi_values[cell_index] <- 1 /
      sum(batch_probabilities^2)
  }
  
  names(ilisi_values) <-
    rownames(knn_result$indices)
  
  ilisi_values
}


calculate_kbet <- function(
    knn_indices,
    batch_labels,
    neighbourhood_size = 50,
    alpha = 0.05,
    number_of_repeats = 100,
    test_fraction = 0.10,
    seed = 2026
) {
  
  if (
    neighbourhood_size - 1 >
    ncol(knn_indices)
  ) {
    stop("kBET neighbourhood exceeds available neighbours.")
  }
  
  batch_factor <- factor(
    batch_labels
  )
  
  batch_codes <- as.integer(
    batch_factor
  )
  
  number_of_batches <- nlevels(
    batch_factor
  )
  
  global_batch_proportions <- as.numeric(
    prop.table(
      table(batch_factor)
    )
  )
  
  expected_counts <-
    neighbourhood_size *
    global_batch_proportions
  
  p_values <- numeric(
    nrow(knn_indices)
  )
  
  for (
    cell_index in
    seq_len(nrow(knn_indices))
  ) {
    
    neighbourhood_cells <- c(
      cell_index,
      knn_indices[
        cell_index,
        seq_len(neighbourhood_size - 1)
      ]
    )
    
    observed_counts <- tabulate(
      batch_codes[neighbourhood_cells],
      nbins = number_of_batches
    )
    
    chi_square_statistic <- sum(
      (
        observed_counts -
          expected_counts
      )^2 /
        expected_counts
    )
    
    p_values[cell_index] <- pchisq(
      chi_square_statistic,
      df = number_of_batches - 1,
      lower.tail = FALSE
    )
  }
  
  rejected <- p_values < alpha
  
  test_size <- ceiling(
    nrow(knn_indices) *
      test_fraction
  )
  
  set.seed(seed)
  
  repeated_rejection_rates <- replicate(
    number_of_repeats,
    {
      tested_cells <- sample.int(
        n = nrow(knn_indices),
        size = test_size,
        replace = FALSE
      )
      
      mean(
        rejected[tested_cells],
        na.rm = TRUE
      )
    }
  )
  
  list(
    summary = data.frame(
      neighbourhood_size =
        neighbourhood_size,
      test_size = test_size,
      number_of_repeats =
        number_of_repeats,
      overall_rejection_rate =
        mean(rejected, na.rm = TRUE),
      mean_repeated_rejection_rate =
        mean(repeated_rejection_rates),
      CI95_lower = unname(
        quantile(
          repeated_rejection_rates,
          0.025
        )
      ),
      CI95_upper = unname(
        quantile(
          repeated_rejection_rates,
          0.975
        )
      )
    ),
    p_values = p_values,
    repeated_rejection_rates =
      repeated_rejection_rates
  )
}


calculate_graph_connectivity <- function(
    knn_indices,
    biological_labels,
    graph_k = 20
) {
  
  if (graph_k > ncol(knn_indices)) {
    stop("graph_k exceeds the available neighbours.")
  }
  
  number_of_cells <- nrow(
    knn_indices
  )
  
  edge_matrix <- cbind(
    rep(
      seq_len(number_of_cells),
      each = graph_k
    ),
    as.vector(
      t(
        knn_indices[
          ,
          seq_len(graph_k),
          drop = FALSE
        ]
      )
    )
  )
  
  knn_graph <- igraph::graph_from_edgelist(
    edge_matrix,
    directed = FALSE
  )
  
  knn_graph <- igraph::simplify(
    knn_graph,
    remove.multiple = TRUE,
    remove.loops = TRUE
  )
  
  biological_levels <- sort(
    unique(
      as.character(biological_labels)
    )
  )
  
  per_population_connectivity <- bind_rows(
    lapply(
      biological_levels,
      function(population_name) {
        
        population_cells <- which(
          biological_labels ==
            population_name
        )
        
        population_graph <-
          igraph::induced_subgraph(
            knn_graph,
            vids = population_cells
          )
        
        component_sizes <-
          igraph::components(
            population_graph
          )$csize
        
        connectivity_value <-
          max(component_sizes) /
          length(population_cells)
        
        data.frame(
          biological_label =
            population_name,
          number_of_cells =
            length(population_cells),
          graph_connectivity =
            connectivity_value,
          stringsAsFactors = FALSE
        )
      }
    )
  )
  
  list(
    overall_mean = mean(
      per_population_connectivity$
        graph_connectivity
    ),
    weighted_mean = weighted.mean(
      per_population_connectivity$
        graph_connectivity,
      per_population_connectivity$
        number_of_cells
    ),
    per_population =
      per_population_connectivity
  )
}

cat(
  "Neighbourhood metric functions defined successfully.\n"
)





# ============================================================
# 8B. Pilot neighbourhood metrics: unintegrated PCA
# ============================================================

cat(
  "Building unintegrated PCA neighbourhood graph...\n"
)

unintegrated_embedding <- getElement(
  workflow_results,
  "Unintegrated PCA"
)$embedding

unintegrated_knn <- build_knn_result(
  embedding_matrix = unintegrated_embedding,
  requested_neighbors = 90
)

cat("Calculating unintegrated iLISI...\n")

unintegrated_ilisi <- calculate_ilisi(
  knn_result = unintegrated_knn,
  batch_labels =
    benchmark_metadata$donor_id_analysis,
  perplexity = 30
)

cat("Calculating unintegrated kBET...\n")

unintegrated_kbet <- calculate_kbet(
  knn_indices =
    unintegrated_knn$indices,
  batch_labels =
    benchmark_metadata$donor_id_analysis,
  neighbourhood_size = 50,
  alpha = 0.05,
  number_of_repeats = 100,
  test_fraction = 0.10,
  seed = benchmark_seed
)

cat(
  "Calculating unintegrated graph connectivity...\n"
)

unintegrated_connectivity <-
  calculate_graph_connectivity(
    knn_indices =
      unintegrated_knn$indices,
    biological_labels =
      benchmark_metadata$biological_label,
    graph_k = 20
  )

number_of_donors <- length(
  unique(
    benchmark_metadata$donor_id_analysis
  )
)

unintegrated_metric_summary <- data.frame(
  workflow = "Unintegrated PCA",
  mean_iLISI = mean(
    unintegrated_ilisi,
    na.rm = TRUE
  ),
  median_iLISI = median(
    unintegrated_ilisi,
    na.rm = TRUE
  ),
  scaled_median_iLISI = (
    median(
      unintegrated_ilisi,
      na.rm = TRUE
    ) - 1
  ) / (
    number_of_donors - 1
  ),
  kBET_rejection_rate =
    unintegrated_kbet$summary$
    mean_repeated_rejection_rate,
  kBET_CI95_lower =
    unintegrated_kbet$summary$
    CI95_lower,
  kBET_CI95_upper =
    unintegrated_kbet$summary$
    CI95_upper,
  graph_connectivity =
    unintegrated_connectivity$
    overall_mean,
  weighted_graph_connectivity =
    unintegrated_connectivity$
    weighted_mean
)

write.csv(
  unintegrated_metric_summary,
  file.path(
    output_dir,
    "pilot_unintegrated_neighbourhood_metrics.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    internal_cell_id = expected_cell_ids,
    iLISI = unintegrated_ilisi,
    kBET_p_value =
      unintegrated_kbet$p_values
  ),
  file.path(
    output_dir,
    "pilot_unintegrated_cell_level_metrics.csv"
  ),
  row.names = FALSE
)

write.csv(
  unintegrated_connectivity$per_population,
  file.path(
    output_dir,
    "pilot_unintegrated_graph_connectivity_by_population.csv"
  ),
  row.names = FALSE
)

print(
  unintegrated_metric_summary %>%
    mutate(
      across(
        where(is.numeric),
        function(x) round(x, 4)
      )
    ),
  row.names = FALSE
)

cat(
  "Unintegrated neighbourhood-metric pilot completed successfully.\n"
)




# ============================================================
# 8C. Define evaluator and calculate CCA neighbourhood metrics
# ============================================================

evaluate_neighbourhood_metrics <- function(
    workflow_name,
    workflow_result
) {
  
  cat(
    "\nEvaluating neighbourhood metrics for",
    workflow_name, "...\n"
  )
  
  knn_result <- build_knn_result(
    embedding_matrix =
      workflow_result$embedding,
    requested_neighbors = 90
  )
  
  ilisi_values <- calculate_ilisi(
    knn_result = knn_result,
    batch_labels =
      benchmark_metadata$donor_id_analysis,
    perplexity = 30
  )
  
  kbet_result <- calculate_kbet(
    knn_indices = knn_result$indices,
    batch_labels =
      benchmark_metadata$donor_id_analysis,
    neighbourhood_size = 50,
    alpha = 0.05,
    number_of_repeats = 100,
    test_fraction = 0.10,
    seed = benchmark_seed
  )
  
  connectivity_result <-
    calculate_graph_connectivity(
      knn_indices = knn_result$indices,
      biological_labels =
        benchmark_metadata$biological_label,
      graph_k = 20
    )
  
  metric_summary <- data.frame(
    workflow = workflow_name,
    mean_iLISI = mean(
      ilisi_values,
      na.rm = TRUE
    ),
    median_iLISI = median(
      ilisi_values,
      na.rm = TRUE
    ),
    scaled_median_iLISI = (
      median(
        ilisi_values,
        na.rm = TRUE
      ) - 1
    ) / (
      number_of_donors - 1
    ),
    kBET_rejection_rate =
      kbet_result$summary$
      mean_repeated_rejection_rate,
    kBET_CI95_lower =
      kbet_result$summary$
      CI95_lower,
    kBET_CI95_upper =
      kbet_result$summary$
      CI95_upper,
    graph_connectivity =
      connectivity_result$overall_mean,
    weighted_graph_connectivity =
      connectivity_result$weighted_mean,
    stringsAsFactors = FALSE
  )
  
  safe_workflow_name <- gsub(
    "[^A-Za-z0-9]+",
    "_",
    workflow_name
  )
  
  write.csv(
    data.frame(
      internal_cell_id =
        expected_cell_ids,
      iLISI = ilisi_values,
      kBET_p_value =
        kbet_result$p_values
    ),
    file.path(
      output_dir,
      paste0(
        safe_workflow_name,
        "_cell_level_neighbourhood_metrics.csv"
      )
    ),
    row.names = FALSE
  )
  
  write.csv(
    data.frame(
      repeat_number =
        seq_along(
          kbet_result$
            repeated_rejection_rates
        ),
      kBET_rejection_rate =
        kbet_result$
        repeated_rejection_rates
    ),
    file.path(
      output_dir,
      paste0(
        safe_workflow_name,
        "_kBET_repeated_rates.csv"
      )
    ),
    row.names = FALSE
  )
  
  connectivity_table <-
    connectivity_result$per_population
  
  connectivity_table$workflow <-
    workflow_name
  
  write.csv(
    connectivity_table,
    file.path(
      output_dir,
      paste0(
        safe_workflow_name,
        "_graph_connectivity_by_population.csv"
      )
    ),
    row.names = FALSE
  )
  
  list(
    summary = metric_summary,
    connectivity =
      connectivity_table
  )
}

cca_neighbourhood_result <-
  evaluate_neighbourhood_metrics(
    workflow_name = "CCA",
    workflow_result =
      getElement(
        workflow_results,
        "CCA"
      )
  )

print(
  cca_neighbourhood_result$summary %>%
    mutate(
      across(
        where(is.numeric),
        function(x) round(x, 4)
      )
    ),
  row.names = FALSE
)

cat("CCA neighbourhood metrics completed and saved.\n")



# ============================================================
# 8D. Calculate RPCA neighbourhood metrics
# ============================================================

rpca_neighbourhood_result <-
  evaluate_neighbourhood_metrics(
    workflow_name = "RPCA",
    workflow_result =
      getElement(
        workflow_results,
        "RPCA"
      )
  )

print(
  rpca_neighbourhood_result$summary %>%
    mutate(
      across(
        where(is.numeric),
        function(x) round(x, 4)
      )
    ),
  row.names = FALSE
)

cat("RPCA neighbourhood metrics completed and saved.\n")




# ============================================================
# 8E. Calculate Harmony neighbourhood metrics
# ============================================================

harmony_neighbourhood_result <-
  evaluate_neighbourhood_metrics(
    workflow_name = "Harmony",
    workflow_result =
      getElement(
        workflow_results,
        "Harmony"
      )
  )

print(
  harmony_neighbourhood_result$summary %>%
    mutate(
      across(
        where(is.numeric),
        function(x) round(x, 4)
      )
    ),
  row.names = FALSE
)

cat("Harmony neighbourhood metrics completed and saved.\n")




# ============================================================
# 9. Assemble the complete expanded benchmark table
# ============================================================

neighbourhood_metrics <- bind_rows(
  unintegrated_metric_summary,
  cca_neighbourhood_result$summary,
  rpca_neighbourhood_result$summary,
  harmony_neighbourhood_result$summary
)

expanded_benchmark_table <- initial_metrics %>%
  left_join(
    neighbourhood_metrics,
    by = "workflow"
  ) %>%
  mutate(
    kBET_acceptance_rate =
      1 - kBET_rejection_rate,
    workflow = factor(
      workflow,
      levels = c(
        "Unintegrated PCA",
        "CCA",
        "RPCA",
        "Harmony"
      )
    )
  ) %>%
  arrange(workflow) %>%
  mutate(
    workflow = as.character(workflow)
  )

stopifnot(
  nrow(expanded_benchmark_table) == 4,
  !anyNA(expanded_benchmark_table)
)

write.csv(
  expanded_benchmark_table,
  file.path(
    output_dir,
    "Table_expanded_integration_benchmark.csv"
  ),
  row.names = FALSE
)

all_connectivity_results <- bind_rows(
  unintegrated_connectivity$per_population %>%
    mutate(
      workflow = "Unintegrated PCA"
    ),
  cca_neighbourhood_result$connectivity,
  rpca_neighbourhood_result$connectivity,
  harmony_neighbourhood_result$connectivity
)

write.csv(
  all_connectivity_results,
  file.path(
    output_dir,
    "Table_graph_connectivity_by_population_and_workflow.csv"
  ),
  row.names = FALSE
)

print(
  expanded_benchmark_table %>%
    select(
      workflow,
      clusters_recovered,
      ARI_vs_full_RPCA_reference,
      NMI_vs_full_RPCA_reference,
      mean_biological_silhouette,
      scaled_median_iLISI,
      kBET_rejection_rate,
      graph_connectivity,
      runtime_seconds
    ) %>%
    mutate(
      across(
        where(is.numeric),
        function(x) round(x, 4)
      )
    ),
  row.names = FALSE
)

cat("Complete expanded benchmark table saved successfully.\n")



# ============================================================
# 10A. Generate final benchmark quality-metric figure
# ============================================================

quality_plot_data <- expanded_benchmark_table %>%
  select(
    workflow,
    ARI_vs_full_RPCA_reference,
    NMI_vs_full_RPCA_reference,
    mean_biological_silhouette,
    scaled_median_iLISI,
    kBET_rejection_rate,
    graph_connectivity
  ) %>%
  pivot_longer(
    cols = -workflow,
    names_to = "metric",
    values_to = "metric_value"
  ) %>%
  mutate(
    workflow = factor(
      workflow,
      levels = c(
        "Unintegrated PCA",
        "CCA",
        "RPCA",
        "Harmony"
      )
    ),
    metric = factor(
      metric,
      levels = c(
        "ARI_vs_full_RPCA_reference",
        "NMI_vs_full_RPCA_reference",
        "mean_biological_silhouette",
        "scaled_median_iLISI",
        "kBET_rejection_rate",
        "graph_connectivity"
      ),
      labels = c(
        "ARI (higher is better)",
        "NMI (higher is better)",
        "Biological ASW (higher is better)",
        "Scaled iLISI (higher is better)",
        "kBET rejection (lower is better)",
        "Graph connectivity (higher is better)"
      )
    )
  )

workflow_colours <- c(
  "Unintegrated PCA" = "#7F7F7F",
  "CCA" = "#377EB8",
  "RPCA" = "#E41A1C",
  "Harmony" = "#4DAF4A"
)

quality_metric_figure <- ggplot(
  quality_plot_data,
  aes(
    x = workflow,
    y = metric_value,
    colour = workflow,
    group = 1
  )
) +
  geom_line(
    colour = "grey75",
    linewidth = 0.6
  ) +
  geom_point(
    size = 4
  ) +
  geom_text(
    aes(
      label = sprintf(
        "%.3f",
        metric_value
      )
    ),
    colour = "black",
    vjust = -1,
    size = 3.6,
    show.legend = FALSE
  ) +
  facet_wrap(
    ~ metric,
    scales = "free_y",
    ncol = 2
  ) +
  scale_colour_manual(
    values = workflow_colours
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0.12, 0.45)
    )
  ) +
  labs(
    title = "Integration workflowbenchmark",
    subtitle = paste0(
      "Common donor-balanced subset: 6,955 cells; ",
      "2,000 HVGs; 20 PCs; k = 20; resolution = 0.6\n",
      "Facet-specific y-axis ranges are used."
    ),
    x = NULL,
    y = "Metric value",
    colour = "Workflow"
  ) +
  theme_bw(
    base_size = 11
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10
    ),
    strip.text = element_text(
      face = "bold",
      size = 11
    ),
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    ),
    legend.position = "bottom",
    panel.grid.minor = element_blank()
  )

print(quality_metric_figure)

ggsave(
  file.path(
    output_dir,
    "Figure_expanded_integration_quality_metrics.png"
  ),
  quality_metric_figure,
  width = 12,
  height = 10,
  dpi = 400,
  bg = "white"
)

ggsave(
  file.path(
    output_dir,
    "Figure_expanded_integration_quality_metrics.pdf"
  ),
  quality_metric_figure,
  width = 12,
  height = 10,
  bg = "white"
)

cat("Final benchmark quality figure saved successfully.\n")


# ============================================================
# 10B. Generate workflow runtime comparison figure
# ============================================================

runtime_plot_data <- expanded_benchmark_table %>%
  select(
    workflow,
    runtime_seconds
  ) %>%
  mutate(
    workflow = factor(
      workflow,
      levels = c(
        "Unintegrated PCA",
        "CCA",
        "RPCA",
        "Harmony"
      )
    )
  )

runtime_figure <- ggplot(
  runtime_plot_data,
  aes(
    x = workflow,
    y = runtime_seconds,
    fill = workflow
  )
) +
  geom_col(
    width = 0.68,
    colour = "black",
    linewidth = 0.4
  ) +
  geom_text(
    aes(
      label = paste0(
        round(runtime_seconds, 1),
        " s"
      )
    ),
    vjust = -0.6,
    size = 4
  ) +
  scale_fill_manual(
    values = workflow_colours
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(0, 0.15)
    )
  ) +
  labs(
    title = "End-to-end workflow runtime",
    subtitle = paste0(
      "Common donor-balanced subset of 6,955 cells; ",
      "identical preprocessing and clustering settings"
    ),
    x = NULL,
    y = "Runtime (seconds)"
  ) +
  theme_bw(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5,
      size = 16
    ),
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 10
    ),
    axis.text.x = element_text(
      angle = 20,
      hjust = 1
    ),
    legend.position = "none",
    panel.grid.minor = element_blank()
  )

print(runtime_figure)

ggsave(
  file.path(
    output_dir,
    "Figure_workflow_runtime_comparison.png"
  ),
  runtime_figure,
  width = 8,
  height = 6,
  dpi = 400,
  bg = "white"
)

ggsave(
  file.path(
    output_dir,
    "Figure_workflow_runtime_comparison.pdf"
  ),
  runtime_figure,
  width = 8,
  height = 6,
  bg = "white"
)

cat("Workflow runtime figure saved successfully.\n")


# ============================================================
# 11. Save metric definitions and analytical summary
# ============================================================

metric_definitions <- data.frame(
  metric = c(
    "ARI",
    "NMI",
    "Biological ASW",
    "Scaled median iLISI",
    "kBET rejection rate",
    "Graph connectivity",
    "Runtime"
  ),
  calculation = c(
    paste0(
      "Agreement between each workflow clustering and the ",
      "retained full-data RPCA resolution-0.6 partition."
    ),
    paste0(
      "Normalized mutual information between each workflow ",
      "clustering and the retained full-data RPCA partition."
    ),
    paste0(
      "Mean silhouette width for marker-supported biological ",
      "labels using 3,000 identically sampled cells."
    ),
    paste0(
      "Entropy-weighted local inverse Simpson index using ",
      "perplexity 30 and 3 x perplexity nearest neighbours; ",
      "median scaled to the interval 0-1 for five donors."
    ),
    paste0(
      "Pearson chi-square donor-composition test using ",
      "50-cell neighbourhoods, alpha 0.05, 10% tested cells, ",
      "and 100 repeated samples."
    ),
    paste0(
      "Mean fraction of cells contained in the largest connected ",
      "component within each marker-supported biological population, ",
      "using a 20-nearest-neighbour graph."
    ),
    paste0(
      "Sequential end-to-end runtime for preprocessing, integration, ",
      "neighbour-graph construction, and clustering."
    )
  ),
  preferred_direction = c(
    "Higher",
    "Higher",
    "Higher",
    "Higher",
    "Lower",
    "Higher",
    "Lower"
  ),
  interpretation_constraint = c(
    "Agreement measure; not independent biological ground truth.",
    "Agreement measure; not independent biological ground truth.",
    "Assesses preservation of broad marker-supported populations.",
    "Assesses local donor mixing.",
    "Assesses residual local donor-composition deviation.",
    "Assesses within-population structural connectivity.",
    "Hardware- and implementation-dependent."
  ),
  stringsAsFactors = FALSE
)

write.csv(
  metric_definitions,
  file.path(
    output_dir,
    "Table_metric_definitions_and_directions.csv"
  ),
  row.names = FALSE
)

summary_file <- file.path(
  output_dir,
  "expanded_integration_benchmark_summary.txt"
)

capture.output(
  {
    cat("EXPANDED INTEGRATION BENCHMARK SUMMARY\n")
    cat("======================================\n\n")
    
    cat("Benchmark design\n")
    cat("----------------\n")
    cat(
      "Four workflows were evaluated on the identical ",
      "donor-balanced subset of 6,955 cells.\n"
    )
    cat(
      "All workflows used 2,000 HVGs, 20 PCs, ",
      "neighbour k = 20, resolution = 0.6, and seed 2026.\n\n"
    )
    
    cat("Complete results\n")
    cat("----------------\n")
    
    print(
      expanded_benchmark_table,
      row.names = FALSE
    )
    
    cat("\nInterpretation\n")
    cat("--------------\n")
    cat(
      "All integrated workflows improved donor mixing relative ",
      "to unintegrated PCA while retaining high graph connectivity.\n"
    )
    cat(
      "RPCA produced the highest scaled median iLISI (0.4490), ",
      "the lowest kBET rejection rate among integrated workflows ",
      "(0.8687), and the highest NMI agreement (0.7231).\n"
    )
    cat(
      "Harmony produced the highest ARI agreement (0.5871), ",
      "the highest graph connectivity (0.9909), and the shortest ",
      "runtime among integrated workflows (10.5 seconds).\n"
    )
    cat(
      "CCA produced the highest biological silhouette value ",
      "(0.4194) but required the longest runtime (82.3 seconds).\n"
    )
    cat(
      "No method dominated every metric. RPCA was retained as the ",
      "primary workflow because it provided the strongest overall ",
      "balance between donor mixing, clustering agreement, biological ",
      "structure preservation, and computational cost.\n\n"
    )
    
    cat("Interpretive limitations\n")
    cat("------------------------\n")
    cat(
      "ARI and NMI quantify agreement with the retained full-data ",
      "RPCA partition and must not be interpreted as accuracy against ",
      "independent biological ground truth.\n"
    )
    cat(
      "The high kBET rejection rates indicate that residual donor ",
      "structure remains after integration and should be reported ",
      "rather than concealed.\n"
    )
    cat(
      "These findings apply to the five healthy synovial donors ",
      "evaluated here and do not establish universal superiority ",
      "across tissues, diseases, platforms, or populations.\n"
    )
  },
  file = summary_file
)

cat("Metric-definition table saved successfully.\n")
cat("Benchmark summary saved:", file.exists(summary_file), "\n")




# ============================================================
# 12. Save session information and verify outputs
# ============================================================

session_information_file <- file.path(
  output_dir,
  "expanded_integration_benchmark_session_info.txt"
)

capture.output(
  sessionInfo(),
  file = session_information_file
)

required_output_files <- c(
  "Table_expanded_integration_benchmark.csv",
  "Table_metric_definitions_and_directions.csv",
  "Table_graph_connectivity_by_population_and_workflow.csv",
  "cell_level_assignments_all_workflows.csv",
  "Figure_expanded_integration_quality_metrics.png",
  "Figure_expanded_integration_quality_metrics.pdf",
  "Figure_workflow_runtime_comparison.png",
  "Figure_workflow_runtime_comparison.pdf",
  "expanded_integration_benchmark_summary.txt",
  "expanded_integration_benchmark_session_info.txt",
  "unintegrated_PCA_result.rds",
  "CCA_integration_result.rds",
  "RPCA_integration_result.rds",
  "Harmony_integration_result.rds"
)

required_output_paths <- file.path(
  output_dir,
  required_output_files
)

missing_required_outputs <- required_output_files[
  !file.exists(required_output_paths)
]

if (length(missing_required_outputs) > 0) {
  stop(
    "Missing required outputs: ",
    paste(
      missing_required_outputs,
      collapse = ", "
    )
  )
}

all_output_paths <- list.files(
  output_dir,
  full.names = TRUE,
  recursive = TRUE
)

output_manifest <- data.frame(
  file_name = basename(all_output_paths),
  relative_path = sub(
    paste0("^", output_dir, "/"),
    "",
    all_output_paths
  ),
  file_size_bytes = file.info(
    all_output_paths
  )$size,
  md5_checksum = unname(
    tools::md5sum(
      all_output_paths
    )
  ),
  stringsAsFactors = FALSE
)

write.csv(
  output_manifest,
  file.path(
    output_dir,
    "expanded_benchmark_output_manifest.csv"
  ),
  row.names = FALSE
)

print(
  output_manifest %>%
    transmute(
      file_name,
      file_size_MB = round(
        file_size_bytes / 1024^2,
        3
      )
    ),
  row.names = FALSE
)

cat(
  "Required outputs verified:",
  length(required_output_files),
  "of",
  length(required_output_files),
  "\n"
)

cat(
  "Session information saved:",
  file.exists(session_information_file),
  "\n"
)

cat(
  "Output manifest saved:",
  file.exists(
    file.path(
      output_dir,
      "expanded_benchmark_output_manifest.csv"
    )
  ),
  "\n"
)