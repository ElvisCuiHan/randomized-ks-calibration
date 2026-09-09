#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  packages <- c("scDesign3", "SingleCellExperiment", "gamlss.dist")
})

cat("R session information\n")
cat("=====================\n")
print(sessionInfo())

cat("\nPackage versions used by scripts/ks_scdesign3_example.R\n")
cat("=======================================================\n")
for (pkg in packages) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("%s: %s\n", pkg, as.character(utils::packageVersion(pkg))))
  } else {
    cat(sprintf("%s: not installed\n", pkg))
  }
}
