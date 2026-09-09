source("scripts/ks_conditional_calibration.R")
suppressPackageStartupMessages(library(SingleCellExperiment))

scenarios <- conditional_scenarios()
for (scenario_index in c(1L, 4L, 3L)) {
  scenario <- scenarios[scenario_index, ]
  design <- conditional_design(scenario)
  set.seed(827)
  counts <- rnbinom(scenario$n, mu = design$mu, size = design$size)
  fast <- conditional_fit(counts, design, scenario$formula)
  sce <- SingleCellExperiment(list(counts = matrix(counts, 1,
         dimnames = list("gene", paste0("c", seq_along(counts))))))
  colData(sce)$pseudotime <- design$t
  data <- scDesign3::construct_data(sce, assay_use = "counts", celltype = NULL,
    pseudotime = "pseudotime", spatial = NULL, other_covariates = NULL, corr_by = "1")
  marginal <- scDesign3::fit_marginal(data, mu_formula = scenario$formula,
    sigma_formula = "1", family_use = "nb", n_cores = 1, usebam = FALSE, simplify = TRUE)
  parameters <- scDesign3::extract_para(sce, marginal_list = marginal, n_cores = 1,
    family_use = "nb", new_covariate = data$new_covariate, data = data$dat)
  stopifnot(isTRUE(all.equal(fast$mu, as.numeric(parameters$mean_mat))),
            isTRUE(all.equal(rep(1 / fast$size, scenario$n), as.numeric(parameters$sigma_mat))))
}
u <- conditional_pit(counts, fast, runif(scenario$n))
statistic <- conditional_statistics(u, design$groups)
expected <- sqrt(length(u)) * unname(ks.test(u, "punif")$statistic)
stopifnot(isTRUE(all.equal(unname(statistic["marginal"]), expected)))
a <- conditional_trial(scenario, 3L, 1L, 3L, 20260908L)
b <- conditional_trial(scenario, 3L, 1L, 5L, 20260908L)
stopifnot(identical(a$observed_u, b$observed_u), identical(a$counts, b$counts),
          identical(a$reference, b$reference[1:3, , ]),
          a$fit_calls == 4L,
          all(a$p == (1 + a$exceedances) / 4))
permutation <- rev(seq_along(u))
stopifnot(isTRUE(all.equal(conditional_statistics(u, design$groups),
  conditional_statistics(u[permutation], design$groups[permutation]))))
separated <- seq(0.001, 0.999, length.out = 240)
groups <- rep(1:4, each = 60)
stopifnot(isTRUE(all.equal(unname(conditional_statistics(separated, groups)["conditional"]), sqrt(30))))
serial <- lapply(1:2, function(index) conditional_trial(scenario, 3L, index, 2L, 20260908L))
parallel <- parallel::mclapply(2:1, function(index) conditional_trial(scenario, 3L, index, 2L, 20260908L),
                             mc.cores = 2, mc.set.seed = FALSE)
stopifnot(identical(serial[[1]]$reference, parallel[[2]]$reference),
          identical(serial[[2]]$reference, parallel[[1]]$reference))
cat("Conditional tests passed: three-model extraction equivalence, KS maximum, seeds, prefixes and worker order\n")
