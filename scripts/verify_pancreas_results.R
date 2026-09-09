suppressPackageStartupMessages(library(stats))
folder <- "results/scdesign3_example"
rows <- read.csv(file.path(folder, "pancreas_refitted_margin_bootstrap.csv"))
stopifnot(nrow(rows) == 10L)
for (i in seq_len(nrow(rows))) {
  row <- rows[i, ]
  draws <- readRDS(file.path(folder, "replicates", paste0(row$gene, ".rds")))
  stopifnot(length(draws$bootstrap_d) == 1999L, all(is.finite(draws$bootstrap_d)))
  k <- sum(draws$bootstrap_d >= draws$observed_d)
  stopifnot(k == row$exceedances)
  stopifnot(isTRUE(all.equal(draws$observed_d, row$observed_ks_D)))
  interval <- binom.test(k, 1999)$conf.int
  stopifnot(isTRUE(all.equal(as.numeric(interval), c(row$tail_mc_lower, row$tail_mc_upper))))
  stopifnot(isTRUE(all.equal(as.numeric(quantile(draws$bootstrap_d, .95, type = 1)), row$bootstrap_q95)))
}
cat("Verified all 19,990 pancreatic bootstrap draws, exceedances, intervals and cutoffs\n")
