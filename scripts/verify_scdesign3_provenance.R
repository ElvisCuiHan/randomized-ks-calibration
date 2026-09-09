suppressPackageStartupMessages({
  library(scDesign3)
  library(SingleCellExperiment)
})
data(example_sce, package = "scDesign3")
path <- "tmp/pdfs/yl_revision_20260907/E-MTAB-3929.sdrf.txt"
if (!file.exists(path)) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  download.file("https://www.ebi.ac.uk/biostudies/files/E-MTAB-3929/E-MTAB-3929.sdrf.txt", path)
}
sdrf <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
source_rows <- match(colnames(example_sce), sdrf[["Source Name"]])
stopifnot(!anyNA(source_rows))
stages <- as.character(colData(example_sce)$cell_type)
stopifnot(identical(stages, sdrf[["Characteristics[developmental stage]"]][source_rows]))
output <- "results/scdesign3_example/embryo_cell_provenance.csv"
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
write.csv(data.frame(
  cell = colnames(example_sce), developmental_stage = stages,
  study = "E-MTAB-3929", ena_sample = sdrf[["Comment[ENA_SAMPLE]"]][source_rows]
), output, row.names = FALSE)
cat("Verified", ncol(example_sce), "cell identifiers and developmental stages against E-MTAB-3929\n")
cat("Bundled object has", nrow(example_sce), "genes and", length(unique(stages)), "embryonic days, not cell types\n")
cat("This identifies the source cells; it does not reconstruct count preprocessing or pseudotime estimation.\n")
