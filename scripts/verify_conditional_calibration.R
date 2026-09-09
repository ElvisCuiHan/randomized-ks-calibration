source("scripts/ks_conditional_calibration.R")

folder <- "results/conditional_calibration"
metadata <- readRDS(file.path(folder, "protocol.rds"))
stopifnot(metadata$mc == 500L, metadata$boot == 199L)
summary <- read.csv(file.path(folder, "summary.csv"))
checked <- list()
total_refits <- 0L
for (i in seq_len(nrow(metadata$scenarios))) {
  scenario <- metadata$scenarios[i, ]
  design <- conditional_design(scenario)
  records <- lapply(seq_len(metadata$mc), function(trial) {
    record <- readRDS(file.path(folder, scenario$scenario, sprintf("trial-%04d.rds", trial)))
    stopifnot(record$boot == metadata$boot, record$n == scenario$n,
              record$formula == scenario$formula, record$failed_fits == 0L,
              record$fit_calls == metadata$boot + 1L,
              all(is.finite(record$reference)), all(dim(record$reference) == c(199, 2, 2)))
    set.seed(record$observed_seed)
    counts <- rnbinom(scenario$n, mu = design$mu, size = design$size)
    u <- conditional_pit(counts, record$fit, runif(scenario$n))
    stopifnot(identical(counts, record$counts), identical(u, record$observed_u),
              isTRUE(all.equal(conditional_statistics(u, design$groups), record$observed)))
    for (statistic in c("marginal", "conditional")) {
      for (reference in c("fixed", "refitted")) {
        k <- sum(record$reference[, statistic, reference] >= record$observed[statistic])
        stopifnot(k == record$exceedances[statistic, reference],
                  (1 + k) / 200 == record$p[statistic, reference])
      }
    }
    record
  })
  expected <- summarize_conditional_trials(records)
  actual <- summary[summary$scenario == scenario$scenario, ]
  rownames(actual) <- NULL
  stopifnot(isTRUE(all.equal(expected, actual, check.attributes = FALSE)))
  total_refits <- total_refits + length(records) * metadata$boot
  checked[[i]] <- expected
}
stopifnot(nrow(summary) == 16L, total_refits == 398000L)
cat("Verified 2000 outer samples, 398000 bootstrap refits and 16 summary rows\n")
