Sys.setenv(SCDESIGN3_FUNCTIONS_ONLY = "1")
source("scripts/ks_scdesign3_example.R")
data(example_count, package = "scDesign3")
data(pseudotime, package = "scDesign3")
genes <- rownames(example_count)[1:2]
counts <- as.matrix(example_count)
start <- proc.time()[["elapsed"]]
a <- run_refitted_margin_bootstrap(counts, pseudotime, genes, boot = 3, workers = 1,
  output_dir = "tmp/pdfs/yl_revision_20260907/r_test_a")
b <- run_refitted_margin_bootstrap(counts, pseudotime, rev(genes), boot = 5, workers = 2,
  output_dir = "tmp/pdfs/yl_revision_20260907/r_test_b")
b <- b[match(a$gene, b$gene), ]
stopifnot(identical(a$observed_ks_D, b$observed_ks_D))
for (gene in genes) {
  x <- readRDS(file.path("tmp/pdfs/yl_revision_20260907/r_test_a/replicates", paste0(gene, ".rds")))
  y <- readRDS(file.path("tmp/pdfs/yl_revision_20260907/r_test_b/replicates", paste0(gene, ".rds")))
  stopifnot(identical(x$observed_u, y$observed_u))
  stopifnot(identical(x$bootstrap_d, y$bootstrap_d[1:3]))
}
cat("Observed PIT and replicate-prefix invariance passed; elapsed seconds:",
    proc.time()[["elapsed"]] - start, "\n")
