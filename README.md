# Multi-Donor scRNA-seq Reproducibility Benchmark
[Download the v1.0.0 benchmark package, saved results, scripts, and companion data files](https://github.com/ersahil2011/multi-donor-scrna-reproducibility-benchmark/releases/tag/v1.0.0).
## Overview

This repository contains R scripts, result tables, figures, and session information supporting the manuscript:

**Computational Benchmarking Framework for Reproducible Multi-Donor Single-Cell RNA-Seq Analysis**

The study evaluates within-dataset clustering stability, donor-representation sensitivity, marker-supported population preservation, parameter sensitivity, and integration-workflow behaviour using five healthy synovial tissue donors.

## Dataset

The analysis uses publicly available scRNA-seq data from GEO accession [GSE246416](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE246416).

The five analysed healthy-control samples are:

- GSM7868582_HC002
- GSM7868583_HC004
- GSM7868584_HCBX003
- GSM7868585_HCBXN001
- GSM7868586_HCBXN005

The study does not perform disease-group comparison, biomarker discovery, or disease-mechanism inference.

## Repository structure

- `scripts/`: final R scripts used for the reported analyses
- `results/`: summary tables, metric outputs, seeds, checksums, and manifests
- `figures/main/`: figures used in the manuscript
- `figures/supplementary/`: supporting figures
- `session_info/`: package and software session information

## Main analyses

1. Reproducibility audit of retained Seurat objects
2. Resolution-wise clustering-stability evaluation
3. Fifty donor-balanced RPCA reruns
4. Canonical-marker-supported annotation validation
5. Matched-size unbalanced sampling control
6. HVG, principal-component, and neighbour-k sensitivity analysis
7. Integration benchmarking of unintegrated PCA, CCA, RPCA, and Harmony
8. Marker-supported population preservation across balanced reruns

## Recommended script order

1. `00_reproducibility_audit.R`
2. `01_quantitative_clustering_stability.R`
3. `03_repeated_donor_balanced_sensitivity.R`
4. `10_marker_annotation_validation.R`
5. `11_matched_size_unbalanced_control.R`
6. `12_parameter_sensitivity.R`
7. `13_expanded_integration_benchmark.R`
8. `14_recalculate_marker_supported_preservation.R`

The scripts `07_generate_figure4_resolution_stability.R` and `09_generate_session_info.R` generate supporting figures and session records.

## Required analytical inputs

The scripts operate on documented intermediate Seurat objects, including:

- `sce_raw_import.rds`
- `sce_qc_final.rds`
- `sce_after_doublet_removal.rds`
- `sce_integrated_layers_rpca.rds`
- `sce_final_clusters_locked_rpca_pcs20_res06.rds`
- `sce_labeled_majorityvote.rds`
- `markers_all_clusters_res06.csv`

Large Seurat objects and raw expression matrices are not stored in this GitHub repository. Object summaries, input checksums where available, analytical settings, seeds, and derived outputs are provided to document the evaluated workflow.

## Principal settings

- Log normalization scale factor: 10,000
- Highly variable genes: 2,000
- Principal components: 20
- Nearest neighbours: 20
- Clustering resolution: 0.6
- Donor-balanced sample size: 1,391 cells per donor; 6,955 cells per iteration
- Donor-balanced repetitions: 50

## Software

The principal revision analyses used R 4.5.1, Seurat 5.3.0, SeuratObject 5.2.0, Harmony 2.0.5, RANN 2.6.2, igraph 2.1.4, mclust 6.1.1, and aricode 1.1.0. Exact session information is provided in `session_info/`.

## Interpretation scope

The retained resolution-0.6 partition is used as an internal reference for within-dataset comparison and is not treated as independent biological ground truth. Results should not be generalized beyond the evaluated healthy synovial tissue dataset without external validation.

## Code availability

This repository is currently maintained privately while institutional, supervisory, and co-author approval for public release is obtained.
