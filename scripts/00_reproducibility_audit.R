# ============================================================
# 00_reproducibility_audit_corrected.R
# Purpose:
# Recover all available methodological and result information
# from saved scRNA-seq RDS objects and CSV output files.
#
# This script only reads existing files and generates summaries.
# It does not modify the original RDS objects.
# ============================================================


# ============================================================
# 1. Confirm required packages are loaded
# ============================================================

required_packages <- c(
  "Seurat",
  "SeuratObject",
  "SingleCellExperiment",
  "dplyr",
  "readr",
  "tibble"
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(
      paste0(
        "Package '", pkg, "' is not available. ",
        "Please install/load it before running this script."
      )
    )
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(SingleCellExperiment)
  library(dplyr)
  library(readr)
  library(tibble)
})

cat("Required packages are available and loaded.\n")


# ============================================================
# 2. Create a fresh output directory
# ============================================================

output_dir <- "audit_outputs"

if (dir.exists(output_dir)) {
  unlink(output_dir, recursive = TRUE)
}

dir.create(output_dir)

cat("Fresh output folder created:", output_dir, "\n")


# ============================================================
# 3. Helper functions
# ============================================================

write_output <- function(file_name, content) {
  writeLines(
    text = content,
    con = file.path(output_dir, file_name)
  )
}

safe_table_text <- function(x, title) {
  c(
    "",
    title,
    paste(rep("-", nchar(title)), collapse = ""),
    capture.output(print(table(x, useNA = "ifany")))
  )
}

get_seurat_assay_names <- function(obj) {
  tryCatch(
    {
      names(obj@assays)
    },
    error = function(e) {
      character(0)
    }
  )
}

get_seurat_reduction_names <- function(obj) {
  tryCatch(
    {
      names(obj@reductions)
    },
    error = function(e) {
      character(0)
    }
  )
}

get_layers_safely <- function(assay_object) {
  tryCatch(
    {
      SeuratObject::Layers(assay_object)
    },
    error = function(e) {
      character(0)
    }
  )
}

get_layer_dimension_safely <- function(obj, assay_name, layer_name) {
  tryCatch(
    {
      mat <- SeuratObject::LayerData(
        object = obj,
        assay = assay_name,
        layer = layer_name
      )
      paste(dim(mat), collapse = " x ")
    },
    error = function(e) {
      paste("Could not extract layer:", conditionMessage(e))
    }
  )
}

get_counts_by_slot_safely <- function(obj, assay_name) {
  tryCatch(
    {
      mat <- Seurat::GetAssayData(
        object = obj,
        assay = assay_name,
        slot = "counts"
      )
      paste(dim(mat), collapse = " x ")
    },
    error = function(e) {
      paste("Could not recover counts through slot access:", conditionMessage(e))
    }
  )
}


# ============================================================
# 4. Record all files present in the working directory
# ============================================================

all_files <- list.files()

write_output(
  "00_available_files.txt",
  c(
    "FILES AVAILABLE IN THE WORKING DIRECTORY",
    "========================================",
    "",
    all_files
  )
)

cat("Available file list recorded.\n")


# ============================================================
# 5. Record R session and package information
# ============================================================

packages_to_check <- c(
  "Seurat",
  "SeuratObject",
  "SingleCellExperiment",
  "SingleR",
  "celldex",
  "clustree",
  "scDblFinder",
  "DoubletFinder",
  "edgeR",
  "limma",
  "muscat",
  "Matrix",
  "dplyr",
  "readr",
  "ggplot2"
)

session_text <- capture.output({
  
  cat("R SESSION INFORMATION\n")
  cat("=====================\n\n")
  print(sessionInfo())
  
  cat("\n\nPACKAGE VERSION CHECK\n")
  cat("=====================\n")
  
  for (pkg in packages_to_check) {
    if (requireNamespace(pkg, quietly = TRUE)) {
      cat(pkg, ":", as.character(packageVersion(pkg)), "\n")
    } else {
      cat(pkg, ": Not installed in the current R environment\n")
    }
  }
})

write_output(
  "01_session_and_package_information.txt",
  session_text
)

cat("Session and package information recorded.\n")


# ============================================================
# 6. Identify all RDS files available
# ============================================================

rds_files_present <- list.files(
  pattern = "\\.rds$",
  ignore.case = TRUE
)

write_output(
  "02_rds_files_detected.txt",
  c(
    "RDS FILES DETECTED IN THE WORKING DIRECTORY",
    "===========================================",
    "",
    if (length(rds_files_present) == 0) {
      "No RDS files found."
    } else {
      rds_files_present
    }
  )
)

cat("RDS files detected:", length(rds_files_present), "\n")


# ============================================================
# 7. Function to summarize each RDS object safely
# ============================================================

summarize_rds_object <- function(file_name) {
  
  cat("Reading object:", file_name, "\n")
  
  base_name <- tools::file_path_sans_ext(basename(file_name))
  summary_file_name <- paste0("object_summary_", base_name, ".txt")
  
  obj <- tryCatch(
    {
      readRDS(file_name)
    },
    error = function(e) {
      write_output(
        summary_file_name,
        c(
          "ERROR WHILE READING OBJECT",
          "==========================",
          "",
          paste("File:", file_name),
          paste("Error:", conditionMessage(e))
        )
      )
      return(NULL)
    }
  )
  
  if (is.null(obj)) {
    return(invisible(NULL))
  }
  
  summary_text <- capture.output({
    
    cat("OBJECT FILE\n")
    cat("===========\n")
    cat(file_name, "\n\n")
    
    cat("OBJECT CLASS\n")
    cat("============\n")
    print(class(obj))
    
    cat("\nOBJECT DIMENSIONS\n")
    cat("=================\n")
    print(dim(obj))
    
    
    # --------------------------------------------------------
    # A. SEURAT OBJECT AUDIT
    # --------------------------------------------------------
    
    if (inherits(obj, "Seurat")) {
      
      cat("\nOBJECT TYPE\n")
      cat("===========\n")
      cat("Seurat object\n")
      
      assay_names <- get_seurat_assay_names(obj)
      reduction_names <- get_seurat_reduction_names(obj)
      
      cat("\nASSAYS AVAILABLE\n")
      cat("================\n")
      print(assay_names)
      
      cat("\nDEFAULT ASSAY\n")
      cat("=============\n")
      print(DefaultAssay(obj))
      
      cat("\nREDUCTIONS AVAILABLE\n")
      cat("====================\n")
      print(reduction_names)
      
      cat("\nNUMBER OF FEATURES AND CELLS\n")
      cat("============================\n")
      cat("Features:", nrow(obj), "\n")
      cat("Cells:", ncol(obj), "\n")
      
      cat("\nMETADATA COLUMNS\n")
      cat("================\n")
      print(colnames(obj@meta.data))
      
      # -----------------------------------
      # Sample and experimental-group fields
      # -----------------------------------
      
      sample_cols <- grep(
        pattern = "orig.ident|sample|group|condition|case|control|status|patient|identity",
        x = colnames(obj@meta.data),
        value = TRUE,
        ignore.case = TRUE
      )
      
      cat("\nPOTENTIAL SAMPLE OR GROUP METADATA COLUMNS\n")
      cat("==========================================\n")
      print(sample_cols)
      
      if (length(sample_cols) > 0) {
        for (col in sample_cols) {
          cat("\nVALUE COUNTS FOR:", col, "\n")
          print(table(obj@meta.data[[col]], useNA = "ifany"))
        }
      }
      
      # -----------------------------------
      # Cluster and resolution fields
      # -----------------------------------
      
      cluster_cols <- grep(
        pattern = "cluster|res\\.|resolution|snn",
        x = colnames(obj@meta.data),
        value = TRUE,
        ignore.case = TRUE
      )
      
      cat("\nPOTENTIAL CLUSTER OR RESOLUTION METADATA COLUMNS\n")
      cat("================================================\n")
      print(cluster_cols)
      
      if (length(cluster_cols) > 0) {
        for (col in cluster_cols) {
          cat("\nVALUE COUNTS FOR:", col, "\n")
          print(table(obj@meta.data[[col]], useNA = "ifany"))
        }
      }
      
      # -----------------------------------
      # Active identity information
      # -----------------------------------
      
      cat("\nACTIVE IDENTITY CLASSES\n")
      cat("=======================\n")
      print(levels(Idents(obj)))
      
      cat("\nACTIVE IDENTITY COUNTS\n")
      cat("======================\n")
      print(table(Idents(obj), useNA = "ifany"))
      
      # -----------------------------------
      # Stored Seurat commands
      # -----------------------------------
      
      cat("\nRECORDED SEURAT COMMANDS\n")
      cat("========================\n")
      
      command_names <- names(obj@commands)
      print(command_names)
      
      if (length(command_names) > 0) {
        for (cmd in command_names) {
          cat("\n---------------------------------------------\n")
          cat("COMMAND:", cmd, "\n")
          print(obj@commands[[cmd]])
        }
      }
      
      # -----------------------------------
      # Assays, layers, and raw counts
      # -----------------------------------
      
      cat("\nASSAY LAYERS AND RAW COUNT AVAILABILITY\n")
      cat("=======================================\n")
      
      if (length(assay_names) == 0) {
        
        cat("No assays were detected in this Seurat object.\n")
        
      } else {
        
        for (assay_name in assay_names) {
          
          cat("\n---------------------------------------------\n")
          cat("ASSAY:", assay_name, "\n")
          
          assay_object <- obj[[assay_name]]
          
          layer_names <- get_layers_safely(assay_object)
          
          cat("Available layers:\n")
          print(layer_names)
          
          if (length(layer_names) > 0) {
            
            count_layers <- layer_names[
              grepl("^counts($|\\.)", layer_names, ignore.case = TRUE)
            ]
            
            cat("Detected counts-related layers:\n")
            print(count_layers)
            
            if (length(count_layers) > 0) {
              
              for (layer_name in count_layers) {
                cat(
                  "Dimension of layer '",
                  layer_name,
                  "': ",
                  get_layer_dimension_safely(obj, assay_name, layer_name),
                  "\n",
                  sep = ""
                )
              }
              
            } else {
              
              cat("No layer beginning with 'counts' was detected.\n")
            }
            
          } else {
            
            cat("Layers() was not available for this assay/version.\n")
            cat("Trying legacy slot-based count recovery:\n")
            cat(get_counts_by_slot_safely(obj, assay_name), "\n")
          }
        }
      }
      
      # -----------------------------------
      # Export metadata
      # -----------------------------------
      
      metadata_export <- obj@meta.data
      metadata_export$cell_barcode <- rownames(metadata_export)
      
      write.csv(
        metadata_export,
        file.path(output_dir, paste0("metadata_", base_name, ".csv")),
        row.names = FALSE
      )
    }
    
    
    # --------------------------------------------------------
    # B. SINGLECELLEXPERIMENT OBJECT AUDIT
    # --------------------------------------------------------
    
    if (inherits(obj, "SingleCellExperiment")) {
      
      cat("\nOBJECT TYPE\n")
      cat("===========\n")
      cat("SingleCellExperiment object\n")
      
      cat("\nNUMBER OF FEATURES AND CELLS\n")
      cat("============================\n")
      cat("Features:", nrow(obj), "\n")
      cat("Cells:", ncol(obj), "\n")
      
      cat("\nASSAYS AVAILABLE\n")
      cat("================\n")
      print(SummarizedExperiment::assayNames(obj))
      
      cat("\nREDUCED DIMENSIONS AVAILABLE\n")
      cat("============================\n")
      print(SingleCellExperiment::reducedDimNames(obj))
      
      sce_metadata <- as.data.frame(SummarizedExperiment::colData(obj))
      
      cat("\nCELL METADATA COLUMNS\n")
      cat("=====================\n")
      print(colnames(sce_metadata))
      
      sample_cols <- grep(
        pattern = "orig.ident|sample|group|condition|case|control|status|patient|identity",
        x = colnames(sce_metadata),
        value = TRUE,
        ignore.case = TRUE
      )
      
      cat("\nPOTENTIAL SAMPLE OR GROUP METADATA COLUMNS\n")
      cat("==========================================\n")
      print(sample_cols)
      
      if (length(sample_cols) > 0) {
        for (col in sample_cols) {
          cat("\nVALUE COUNTS FOR:", col, "\n")
          print(table(sce_metadata[[col]], useNA = "ifany"))
        }
      }
      
      cluster_cols <- grep(
        pattern = "cluster|res\\.|resolution|snn|label|celltype|cell_type",
        x = colnames(sce_metadata),
        value = TRUE,
        ignore.case = TRUE
      )
      
      cat("\nPOTENTIAL CLUSTER OR LABEL METADATA COLUMNS\n")
      cat("===========================================\n")
      print(cluster_cols)
      
      if (length(cluster_cols) > 0) {
        for (col in cluster_cols) {
          cat("\nVALUE COUNTS FOR:", col, "\n")
          print(table(sce_metadata[[col]], useNA = "ifany"))
        }
      }
      
      cat("\nRAW COUNT AVAILABILITY\n")
      cat("======================\n")
      
      sce_assays <- SummarizedExperiment::assayNames(obj)
      
      if ("counts" %in% sce_assays) {
        count_dims <- dim(SummarizedExperiment::assay(obj, "counts"))
        cat(
          "Counts assay available. Dimension:",
          paste(count_dims, collapse = " x "),
          "\n"
        )
      } else {
        cat("A counts assay named 'counts' was not detected.\n")
      }
      
      sce_metadata$cell_barcode <- colnames(obj)
      
      write.csv(
        sce_metadata,
        file.path(output_dir, paste0("metadata_", base_name, ".csv")),
        row.names = FALSE
      )
    }
    
    
    # --------------------------------------------------------
    # C. UNKNOWN OBJECT TYPE
    # --------------------------------------------------------
    
    if (!inherits(obj, "Seurat") &&
        !inherits(obj, "SingleCellExperiment")) {
      
      cat("\nOBJECT TYPE STATUS\n")
      cat("==================\n")
      cat("This object is neither a Seurat object nor a SingleCellExperiment object.\n")
      cat("Further inspection may be required.\n")
    }
    
  })
  
  write_output(summary_file_name, summary_text)
  
  rm(obj)
  invisible(gc())
}


# ============================================================
# 8. Audit every detected RDS file
# ============================================================

if (length(rds_files_present) == 0) {
  
  write_output(
    "03_no_rds_files_found.txt",
    c(
      "NO RDS FILES FOUND",
      "==================",
      "",
      "No .rds file was detected in the current working directory."
    )
  )
  
} else {
  
  for (rds_file in rds_files_present) {
    summarize_rds_object(rds_file)
  }
}

cat("RDS object audit completed.\n")


# ============================================================
# 9. Audit the final cell-label CSV file
# ============================================================

label_files <- list.files(
  pattern = "^cell_labels_final.*\\.csv$",
  ignore.case = TRUE
)

if (length(label_files) > 0) {
  
  label_file <- label_files[1]
  
  labels <- read.csv(
    label_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  label_summary <- capture.output({
    
    cat("CELL LABEL FILE USED\n")
    cat("====================\n")
    cat(label_file, "\n\n")
    
    cat("DIMENSIONS\n")
    cat("==========\n")
    print(dim(labels))
    
    cat("\nCOLUMN NAMES\n")
    cat("============\n")
    print(colnames(labels))
    
    cat("\nFIRST SIX ROWS\n")
    cat("==============\n")
    print(head(labels))
    
    possible_cluster_cols <- grep(
      pattern = "cluster|seurat",
      x = colnames(labels),
      value = TRUE,
      ignore.case = TRUE
    )
    
    possible_label_cols <- grep(
      pattern = "label|celltype|cell_type|annotation|majority|final",
      x = colnames(labels),
      value = TRUE,
      ignore.case = TRUE
    )
    
    cat("\nPOSSIBLE CLUSTER COLUMNS\n")
    cat("========================\n")
    print(possible_cluster_cols)
    
    if (length(possible_cluster_cols) > 0) {
      for (col in possible_cluster_cols) {
        cat("\nCOUNT TABLE FOR CLUSTER COLUMN:", col, "\n")
        print(table(labels[[col]], useNA = "ifany"))
      }
    }
    
    cat("\nPOSSIBLE CELL-LABEL COLUMNS\n")
    cat("===========================\n")
    print(possible_label_cols)
    
    if (length(possible_label_cols) > 0) {
      for (col in possible_label_cols) {
        cat("\nCOUNT TABLE FOR CELL-LABEL COLUMN:", col, "\n")
        print(table(labels[[col]], useNA = "ifany"))
      }
    }
    
    if (length(possible_cluster_cols) > 0 &&
        length(possible_label_cols) > 0) {
      
      cat("\nCROSS-TABULATION OF CLUSTER AND FINAL LABEL\n")
      cat("==========================================\n")
      
      print(
        table(
          labels[[possible_cluster_cols[1]]],
          labels[[possible_label_cols[1]]],
          useNA = "ifany"
        )
      )
    }
    
  })
  
  write_output(
    "04_cell_label_summary.txt",
    label_summary
  )
  
  write.csv(
    labels,
    file.path(output_dir, "cell_labels_checked_copy.csv"),
    row.names = FALSE
  )
  
  cat("Cell-label CSV audit completed.\n")
  
} else {
  
  write_output(
    "04_cell_label_summary.txt",
    c(
      "CELL LABEL FILE NOT FOUND",
      "=========================",
      "",
      "No file matching cell_labels_final*.csv was found."
    )
  )
}


# ============================================================
# 10. Audit marker-gene result file
# ============================================================

marker_files <- list.files(
  pattern = "^markers_all_clusters_res06.*\\.csv$",
  ignore.case = TRUE
)

if (length(marker_files) > 0) {
  
  marker_file <- marker_files[1]
  
  markers <- read.csv(
    marker_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  marker_summary <- capture.output({
    
    cat("MARKER FILE USED\n")
    cat("================\n")
    cat(marker_file, "\n\n")
    
    cat("DIMENSIONS\n")
    cat("==========\n")
    print(dim(markers))
    
    cat("\nCOLUMN NAMES\n")
    cat("============\n")
    print(colnames(markers))
    
    cat("\nFIRST SIX ROWS\n")
    cat("==============\n")
    print(head(markers))
    
    possible_cluster_cols <- grep(
      pattern = "cluster",
      x = colnames(markers),
      value = TRUE,
      ignore.case = TRUE
    )
    
    possible_gene_cols <- grep(
      pattern = "gene|feature",
      x = colnames(markers),
      value = TRUE,
      ignore.case = TRUE
    )
    
    cat("\nPOSSIBLE CLUSTER COLUMNS\n")
    cat("========================\n")
    print(possible_cluster_cols)
    
    cat("\nPOSSIBLE GENE COLUMNS\n")
    cat("=====================\n")
    print(possible_gene_cols)
    
    if (length(possible_cluster_cols) > 0) {
      for (col in possible_cluster_cols) {
        cat("\nNUMBER OF MARKER ROWS PER CLUSTER USING:", col, "\n")
        print(table(markers[[col]], useNA = "ifany"))
      }
    }
    
  })
  
  write_output(
    "05_marker_gene_summary.txt",
    marker_summary
  )
  
  cat("Marker-gene CSV audit completed.\n")
  
} else {
  
  write_output(
    "05_marker_gene_summary.txt",
    c(
      "MARKER FILE NOT FOUND",
      "=====================",
      "",
      "No file matching markers_all_clusters_res06*.csv was found."
    )
  )
}


# ============================================================
# 11. Audit all pseudobulk differential-expression results
# ============================================================

de_all_files <- list.files(
  pattern = "^pseudobulk_DE_all_celltypes_tidy.*\\.csv$",
  ignore.case = TRUE
)

if (length(de_all_files) > 0) {
  
  de_all_file <- de_all_files[1]
  
  de_all <- read.csv(
    de_all_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  de_all_summary <- capture.output({
    
    cat("PSEUDOBULK ALL-RESULT FILE USED\n")
    cat("===============================\n")
    cat(de_all_file, "\n\n")
    
    cat("DIMENSIONS\n")
    cat("==========\n")
    print(dim(de_all))
    
    cat("\nCOLUMN NAMES\n")
    cat("============\n")
    print(colnames(de_all))
    
    cat("\nFIRST SIX ROWS\n")
    cat("==============\n")
    print(head(de_all))
    
    possible_celltype_cols <- grep(
      pattern = "cell|type|label",
      x = colnames(de_all),
      value = TRUE,
      ignore.case = TRUE
    )
    
    possible_gene_cols <- grep(
      pattern = "gene|feature",
      x = colnames(de_all),
      value = TRUE,
      ignore.case = TRUE
    )
    
    possible_logfc_cols <- grep(
      pattern = "logFC|log2FoldChange|logfold|lfc",
      x = colnames(de_all),
      value = TRUE,
      ignore.case = TRUE
    )
    
    possible_fdr_cols <- grep(
      pattern = "fdr|adj|padj",
      x = colnames(de_all),
      value = TRUE,
      ignore.case = TRUE
    )
    
    cat("\nPOSSIBLE CELL-TYPE COLUMNS\n")
    cat("==========================\n")
    print(possible_celltype_cols)
    
    cat("\nPOSSIBLE GENE COLUMNS\n")
    cat("=====================\n")
    print(possible_gene_cols)
    
    cat("\nPOSSIBLE LOGFC COLUMNS\n")
    cat("======================\n")
    print(possible_logfc_cols)
    
    cat("\nPOSSIBLE FDR OR ADJUSTED P-VALUE COLUMNS\n")
    cat("========================================\n")
    print(possible_fdr_cols)
    
    if (length(possible_celltype_cols) > 0) {
      for (col in possible_celltype_cols) {
        cat("\nNUMBER OF TESTED ROWS BY:", col, "\n")
        print(table(de_all[[col]], useNA = "ifany"))
      }
    }
    
  })
  
  write_output(
    "06_pseudobulk_all_results_summary.txt",
    de_all_summary
  )
  
  cat("Complete pseudobulk result audit completed.\n")
  
} else {
  
  write_output(
    "06_pseudobulk_all_results_summary.txt",
    c(
      "PSEUDOBULK ALL-RESULT FILE NOT FOUND",
      "====================================",
      "",
      "No file matching pseudobulk_DE_all_celltypes_tidy*.csv was found."
    )
  )
}


# ============================================================
# 12. Audit significant pseudobulk result file
# ============================================================

de_sig_files <- list.files(
  pattern = "^pseudobulk_DE_significant_by_celltype.*\\.csv$",
  ignore.case = TRUE
)

if (length(de_sig_files) > 0) {
  
  de_sig_file <- de_sig_files[1]
  
  de_sig <- read.csv(
    de_sig_file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  de_sig_summary <- capture.output({
    
    cat("PSEUDOBULK SIGNIFICANT-RESULT FILE USED\n")
    cat("=======================================\n")
    cat(de_sig_file, "\n\n")
    
    cat("DIMENSIONS\n")
    cat("==========\n")
    print(dim(de_sig))
    
    cat("\nCOLUMN NAMES\n")
    cat("============\n")
    print(colnames(de_sig))
    
    cat("\nALL ROWS IN SIGNIFICANT RESULT FILE\n")
    cat("===================================\n")
    print(de_sig)
    
  })
  
  write_output(
    "07_pseudobulk_significant_results_summary.txt",
    de_sig_summary
  )
  
  cat("Significant pseudobulk result audit completed.\n")
  
} else {
  
  write_output(
    "07_pseudobulk_significant_results_summary.txt",
    c(
      "PSEUDOBULK SIGNIFICANT-RESULT FILE NOT FOUND",
      "============================================",
      "",
      "No file matching pseudobulk_DE_significant_by_celltype*.csv was found."
    )
  )
}


# ============================================================
# 13. Create README and completion message
# ============================================================

readme_text <- c(
  "REPRODUCIBILITY AUDIT COMPLETED SUCCESSFULLY",
  "============================================",
  "",
  "This folder contains extracted information from the available",
  "RDS analysis objects and CSV result files.",
  "",
  "Please compress the complete 'audit_outputs' folder into a ZIP file",
  "and upload it for manuscript and benchmarking assessment.",
  "",
  "Most important files include:",
  "01_session_and_package_information.txt",
  "02_rds_files_detected.txt",
  "object_summary_<object name>.txt",
  "metadata_<object name>.csv",
  "04_cell_label_summary.txt",
  "05_marker_gene_summary.txt",
  "06_pseudobulk_all_results_summary.txt",
  "07_pseudobulk_significant_results_summary.txt"
)

write_output(
  "AUDIT_COMPLETE_README.txt",
  readme_text
)

cat("\n====================================================\n")
cat("REPRODUCIBILITY AUDIT COMPLETED SUCCESSFULLY\n")
cat("====================================================\n")
cat("A folder named 'audit_outputs' has been created.\n")
cat("Compress this folder into a ZIP file and upload it here.\n")