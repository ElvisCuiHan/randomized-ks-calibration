#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(scDesign3))

conditional_scenarios <- function() {
  data.frame(
    scenario = c("linear_null", "periodic_misspecified", "periodic_null", "spline_null"),
    n = c(240L, 240L, 240L, 2087L),
    formula = c("pseudotime", "pseudotime",
                "pseudotime + sin(2*pi*pseudotime) + cos(2*pi*pseudotime)",
                "s(pseudotime, bs='cr', k=5)"),
    stringsAsFactors = FALSE
  )
}

conditional_design <- function(scenario) {
  t <- (seq_len(scenario$n) - 0.5) / scenario$n
  eta <- 1 + 1.2 * t
  if (grepl("periodic", scenario$scenario)) eta <- eta + 0.8 * sin(2 * pi * t)
  if (scenario$scenario == "spline_null") {
    basis <- mgcv::smoothCon(mgcv::s(pseudotime, bs = "cr", k = 5),
                             data = data.frame(pseudotime = t))[[1]]$X
    eta <- eta + as.numeric(basis %*% qr.coef(qr(basis), 0.8 * sin(2 * pi * t)))
  }
  list(t = t, mu = exp(eta), size = 3,
       groups = cut(t, quantile(t, seq(0, 1, 0.25)), include.lowest = TRUE,
                    labels = FALSE),
       data = list(count_mat = matrix(0, length(t), 1, dimnames = list(NULL, "gene")),
                   dat = data.frame(pseudotime = t, corr_group = 1), filtered_gene = NULL))
}

conditional_fit <- function(counts, design, formula) {
  data <- design$data
  data$count_mat[, 1] <- counts
  result <- scDesign3::fit_marginal(
    data, mu_formula = formula, sigma_formula = "1", family_use = "nb",
    n_cores = 1L, usebam = FALSE, simplify = TRUE, trace = TRUE
  )$gene
  fit <- result$fit
  if (!inherits(fit, "gam") || !isTRUE(fit$converged) ||
      (!is.null(fit$outer.info) && fit$outer.info$conv != "full convergence")) {
    stop("scDesign3 NB fit did not converge")
  }
  mu <- as.numeric(predict(fit, type = "response"))
  size <- fit$family$getTheta(TRUE)
  if (length(size) != 1L || !is.finite(size) || size <= 0 ||
      any(!is.finite(mu)) || any(mu <= 0)) stop("Invalid fitted NB parameters")
  list(mu = mu, size = size, coefficients = coef(fit), edf = sum(fit$edf),
       smoothing = fit$sp, warnings = result$warning)
}

conditional_pit <- function(counts, fit, randomizers) {
  pnbinom(counts - 1, mu = fit$mu, size = fit$size) +
    randomizers * dnbinom(counts, mu = fit$mu, size = fit$size)
}

conditional_statistics <- function(u, groups) {
  n <- length(u)
  ordered <- sort(u)
  d <- max(seq_len(n) / n - ordered, ordered - (seq_len(n) - 1) / n)
  samples <- split(u, groups)
  pairs <- combn(seq_along(samples), 2)
  pair_statistics <- apply(pairs, 2, function(pair) {
    a <- samples[[pair[1]]]
    b <- samples[[pair[2]]]
    distance <- suppressWarnings(unname(stats::ks.test(a, b, exact = FALSE)$statistic))
    sqrt(length(a) * length(b) / (length(a) + length(b))) * distance
  })
  c(marginal = sqrt(n) * d, conditional = max(pair_statistics))
}

conditional_trial <- function(scenario, scenario_id, trial, boot, seed, output_dir = NULL) {
  observed_seed <- as.integer(seed + scenario_id * 100000L + trial)
  bootstrap_seed <- as.integer(observed_seed + 10000000L)
  path <- if (!is.null(output_dir)) file.path(output_dir, scenario$scenario,
                                             sprintf("trial-%04d.rds", trial)) else NULL
  if (!is.null(path) && file.exists(path)) {
    saved <- readRDS(path)
    stopifnot(saved$boot == boot, saved$observed_seed == observed_seed,
              saved$n == scenario$n, saved$formula == scenario$formula,
              saved$protocol == "20260908-v1")
    return(saved)
  }
  start <- proc.time()[["elapsed"]]
  design <- conditional_design(scenario)
  set.seed(observed_seed)
  counts <- rnbinom(scenario$n, mu = design$mu, size = design$size)
  fit <- conditional_fit(counts, design, scenario$formula)
  u <- conditional_pit(counts, fit, runif(scenario$n))
  observed <- conditional_statistics(u, design$groups)
  reference <- array(NA_real_, c(boot, 2L, 2L),
                     dimnames = list(NULL, names(observed), c("fixed", "refitted")))
  warning_count <- length(fit$warnings)
  set.seed(bootstrap_seed)
  for (b in seq_len(boot)) {
    xb <- rnbinom(scenario$n, mu = fit$mu, size = fit$size)
    fb <- conditional_fit(xb, design, scenario$formula)
    vb <- runif(scenario$n)
    reference[b, , "fixed"] <- conditional_statistics(conditional_pit(xb, fit, vb), design$groups)
    reference[b, , "refitted"] <- conditional_statistics(conditional_pit(xb, fb, vb), design$groups)
    warning_count <- warning_count + length(fb$warnings)
  }
  stopifnot(all(is.finite(reference)))
  exceedances <- apply(reference, c(2, 3), function(values) 0L)
  for (statistic in names(observed)) {
    for (method in c("fixed", "refitted")) {
      exceedances[statistic, method] <- sum(reference[, statistic, method] >= observed[statistic])
    }
  }
  result <- list(protocol = "20260908-v1", scenario = scenario$scenario, trial = trial,
                 n = scenario$n, formula = scenario$formula, boot = boot,
                 observed_seed = observed_seed, bootstrap_seed = bootstrap_seed,
                 counts = counts, observed_u = u, observed = observed, fit = fit,
                 reference = reference, exceedances = exceedances,
                 p = (1 + exceedances) / (boot + 1),
                 failed_fits = 0L, fit_calls = boot + 1L, warning_count = warning_count,
                 elapsed_seconds = proc.time()[["elapsed"]] - start)
  if (!is.null(path)) saveRDS(result, path)
  if (trial %% 50 == 0) message(sprintf("%s: outer sample %d complete", scenario$scenario, trial))
  result
}

summarize_conditional_trials <- function(records) {
  first <- records[[1]]
  output <- list()
  for (statistic in c("marginal", "conditional")) {
    for (method in c("fixed", "refitted")) {
      rejected <- vapply(records, function(x) x$p[statistic, method] <= 0.05, logical(1))
      count <- sum(rejected)
      interval <- binom.test(count, length(records))$conf.int
      output[[length(output) + 1L]] <- data.frame(
        scenario = first$scenario, n = first$n, statistic = statistic, reference = method,
        mc = length(records), boot = first$boot, rejections = count, rate = mean(rejected),
        lower = interval[1], upper = interval[2],
        fit_calls = sum(vapply(records, `[[`, numeric(1), "fit_calls")),
        failed_fits = sum(vapply(records, `[[`, numeric(1), "failed_fits")),
        warning_count = sum(vapply(records, `[[`, numeric(1), "warning_count")),
        summed_worker_seconds = sum(vapply(records, `[[`, numeric(1), "elapsed_seconds"))
      )
    }
  }
  do.call(rbind, output)
}

main <- function() {
  mc <- as.integer(Sys.getenv("CONDITIONAL_MC", "500"))
  boot <- as.integer(Sys.getenv("CONDITIONAL_B", "199"))
  workers <- as.integer(Sys.getenv("CONDITIONAL_WORKERS", "4"))
  seed <- 20260908L
  output_dir <- Sys.getenv("CONDITIONAL_OUTPUT", "results/conditional_calibration")
  stopifnot(mc > 0, boot > 0, workers > 0)
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  scenarios <- conditional_scenarios()
  saveRDS(list(scenarios = scenarios, mc = mc, boot = boot, seed = seed,
               workers = workers, session = sessionInfo()), file.path(output_dir, "protocol.rds"))
  summaries <- list()
  start <- proc.time()[["elapsed"]]
  for (i in seq_len(nrow(scenarios))) {
    scenario <- scenarios[i, ]
    dir.create(file.path(output_dir, scenario$scenario), showWarnings = FALSE)
    records <- parallel::mclapply(seq_len(mc), function(trial) {
      conditional_trial(scenario, i, trial, boot, seed, output_dir)
    }, mc.cores = workers, mc.set.seed = FALSE)
    if (any(vapply(records, inherits, logical(1), "try-error"))) stop("A conditional trial failed")
    summaries[[i]] <- summarize_conditional_trials(records)
    write.csv(do.call(rbind, summaries), file.path(output_dir, "summary.csv"), row.names = FALSE)
    print(summaries[[i]])
  }
  writeLines(c(sprintf("elapsed_wall_seconds=%.3f", proc.time()[["elapsed"]] - start),
               capture.output(sessionInfo())), file.path(output_dir, "session_info.txt"))
}

if (sys.nframe() == 0L) main()
