# ============================================================
# 09_generate_session_info.R
#
# Purpose:
# Generate session_info.txt for anonymous code submission.
# ============================================================

output_file <- "session_info.txt"

sink(output_file)

cat("Session information for reproducibility package\n")
cat("================================================\n\n")

cat("Generated on:\n")
print(Sys.time())

cat("\nR session information:\n")
print(sessionInfo())

cat("\nInstalled package versions relevant to the workflow:\n")

packages_to_report <- c(
  "Seurat",
  "SeuratObject",
  "Matrix",
  "dplyr",
  "ggplot2",
  "tidyr",
  "patchwork",
  "mclust",
  "aricode",
  "cluster",
  "SingleR",
  "celldex",
  "SingleCellExperiment",
  "edgeR",
  "limma",
  "scDblFinder",
  "clustree"
)

for (pkg in packages_to_report) {
  cat("\nPackage:", pkg, "\n")
  if (requireNamespace(pkg, quietly = TRUE)) {
    print(packageVersion(pkg))
  } else {
    cat("Not installed or not available in this environment\n")
  }
}

sink()

cat("session_info.txt generated successfully.\n")