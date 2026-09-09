#!/usr/bin/env Rscript

# Single-cell residual diagnostics using objects distributed with scDesign3.
# The script includes a formal refitted marginal bootstrap for the pancreatic
# endocrinogenesis data. It diagnoses conditional margins and does not benchmark
# the full scDesign3 simulator or fitted copula.

suppressPackageStartupMessages({
  library(scDesign3)
  library(SingleCellExperiment)
  library(gamlss.dist)
})

set.seed(20260609)

fit_nb_mom <- function(x) {
  mu <- max(mean(x), 1e-8)
  vv <- var(as.numeric(x))
  size <- if (is.na(vv) || vv <= mu * 1.02) 1000 else max(0.05, mu^2 / (vv - mu))
  c(size = size, mean = mu)
}

randomized_nb_pit <- function(x, size, mu) {
  left <- ifelse(x == 0, 0, pnbinom(x - 1, size = size, mu = mu))
  mass <- dnbinom(x, size = size, mu = mu)
  pmin(pmax(left + runif(length(x)) * mass, 1e-12), 1 - 1e-12)
}

randomized_nbi_pit <- function(x, mu, sigma) {
  left <- ifelse(x == 0, 0, gamlss.dist::pNBI(x - 1, mu = mu, sigma = sigma))
  mass <- gamlss.dist::dNBI(x, mu = mu, sigma = sigma)
  pmin(pmax(left + runif(length(x)) * mass, 1e-12), 1 - 1e-12)
}

ks_stat_uniform <- function(u) {
  u <- sort(as.numeric(u))
  n <- length(u)
  max(max(seq_len(n) / n - u), max(u - (seq_len(n) - 1) / n))
}

pairwise_ks_stat <- function(x, y) {
  suppressWarnings(as.numeric(stats::ks.test(x, y, exact = FALSE)$statistic))
}

make_u_mat <- function(counts) {
  matrix(
    NA_real_,
    nrow = nrow(counts),
    ncol = ncol(counts),
    dimnames = list(rownames(counts), colnames(counts))
  )
}

summarize_residuals <- function(u_mat, pseudotime) {
  gene_ks <- apply(u_mat, 1, ks_stat_uniform)
  naive_p <- vapply(seq_len(nrow(u_mat)), function(g) {
    suppressWarnings(stats::ks.test(u_mat[g, ], "punif")$p.value)
  }, numeric(1))

  resid_cor <- stats::cor(t(u_mat), method = "pearson")
  max_abs_corr <- max(abs(resid_cor[upper.tri(resid_cor)]))

  breaks <- unique(stats::quantile(pseudotime, probs = seq(0, 1, 0.25), na.rm = TRUE))
  bins <- cut(pseudotime, breaks = breaks, include.lowest = TRUE, labels = FALSE)
  quartile_values <- lapply(sort(unique(bins)), function(b) as.numeric(u_mat[, bins == b]))
  max_quartile_ks <- 0
  for (i in seq_along(quartile_values)) {
    for (j in seq_len(i - 1)) {
      max_quartile_ks <- max(
        max_quartile_ks,
        pairwise_ks_stat(quartile_values[[i]], quartile_values[[j]])
      )
    }
  }

  c(
    median_gene_ks_D = median(gene_ks),
    max_gene_ks_D = max(gene_ks),
    min_naive_uniform_p = min(naive_p),
    max_abs_pairwise_residual_corr = max_abs_corr,
    max_pseudotime_quartile_ks = max_quartile_ks
  )
}

run_diagnostic <- function(dataset_name, sce, celltype_arg = NULL) {
  counts <- as.matrix(assay(sce, "counts"))
  cell_info <- as.data.frame(colData(sce))
  pseudotime <- cell_info$pseudotime

  fits <- t(apply(counts, 1, fit_nb_mom))
  moment_u <- make_u_mat(counts)
  for (g in seq_len(nrow(counts))) {
    moment_u[g, ] <- randomized_nb_pit(counts[g, ], fits[g, "size"], fits[g, "mean"])
  }

  scdesign3_data <- construct_data(
    sce = sce,
    assay_use = "counts",
    celltype = celltype_arg,
    pseudotime = "pseudotime",
    spatial = NULL,
    other_covariates = NULL,
    corr_by = "1"
  )

  scdesign3_marginal <- fit_marginal(
    data = scdesign3_data,
    mu_formula = "s(pseudotime, bs = 'cr', k = 5)",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 1,
    usebam = FALSE,
    simplify = TRUE
  )

  scdesign3_para <- extract_para(
    sce = sce,
    marginal_list = scdesign3_marginal,
    n_cores = 1,
    family_use = "nb",
    new_covariate = scdesign3_data$new_covariate,
    data = scdesign3_data$dat
  )

  mean_mat <- t(as.matrix(scdesign3_para$mean_mat))
  sigma_mat <- t(as.matrix(scdesign3_para$sigma_mat))
  scdesign3_u <- make_u_mat(counts)
  for (g in seq_len(nrow(counts))) {
    scdesign3_u[g, ] <- randomized_nbi_pit(counts[g, ], mean_mat[g, ], sigma_mat[g, ])
  }

  summary_mat <- rbind(
    moment_nb = summarize_residuals(moment_u, pseudotime),
    scdesign3_conditional_nb = summarize_residuals(scdesign3_u, pseudotime)
  )

  cell_types <- if (is.null(celltype_arg)) {
    "NA"
  } else {
    as.character(length(unique(cell_info[[celltype_arg]])))
  }

  list(
    metadata = data.frame(
      dataset = dataset_name,
      cells = ncol(counts),
      genes = nrow(counts),
      cell_types = cell_types,
      pseudotime_min = min(pseudotime, na.rm = TRUE),
      pseudotime_max = max(pseudotime, na.rm = TRUE),
      stringsAsFactors = FALSE
    ),
    summary = data.frame(
      dataset = dataset_name,
      fit = rownames(summary_mat),
      summary_mat,
      row.names = NULL,
      check.names = FALSE
    )
  )
}

fit_single_scdesign3_margin <- function(counts, gene_name, pseudotime) {
  one_sce <- SingleCellExperiment(
    list(counts = matrix(
      as.integer(counts),
      nrow = 1,
      dimnames = list(gene_name, paste0("cell", seq_along(counts)))
    ))
  )
  colData(one_sce)$pseudotime <- pseudotime

  one_data <- construct_data(
    sce = one_sce,
    assay_use = "counts",
    celltype = NULL,
    pseudotime = "pseudotime",
    spatial = NULL,
    other_covariates = NULL,
    corr_by = "1"
  )
  one_marginal <- fit_marginal(
    data = one_data,
    mu_formula = "s(pseudotime, bs = 'cr', k = 5)",
    sigma_formula = "1",
    family_use = "nb",
    n_cores = 1,
    usebam = FALSE,
    simplify = TRUE
  )
  one_para <- extract_para(
    sce = one_sce,
    marginal_list = one_marginal,
    n_cores = 1,
    family_use = "nb",
    new_covariate = one_data$new_covariate,
    data = one_data$dat
  )

  list(
    mu = as.numeric(one_para$mean_mat),
    sigma = as.numeric(one_para$sigma_mat)
  )
}

run_refitted_margin_bootstrap <- function(
  counts,
  pseudotime,
  genes,
  boot = 99,
  seed = 20260815,
  workers = 1L,
  output_dir = "results/scdesign3_example"
) {
  stopifnot(boot > 0L, workers > 0L)
  dir.create(file.path(output_dir, "replicates"), recursive = TRUE, showWarnings = FALSE)
  run_gene <- function(gene_index) {
    gene <- genes[[gene_index]]
    gene_seed <- seed + 1009L * (match(gene, rownames(counts)) - 1L)
    set.seed(gene_seed)
    observed_counts <- as.integer(counts[gene, ])
    start <- proc.time()[["elapsed"]]
    observed_fit <- fit_single_scdesign3_margin(observed_counts, gene, pseudotime)
    observed_u <- randomized_nbi_pit(
      observed_counts,
      observed_fit$mu,
      observed_fit$sigma
    )
    observed_d <- ks_stat_uniform(observed_u)
    naive_p <- suppressWarnings(stats::ks.test(observed_u, "punif")$p.value)

    set.seed(gene_seed + 100000L)
    bootstrap_d <- numeric(boot)
    failed <- 0L
    for (b in seq_len(boot)) {
      bootstrap_counts <- gamlss.dist::rNBI(
        length(observed_counts),
        mu = observed_fit$mu,
        sigma = observed_fit$sigma
      )
      bootstrap_fit <- tryCatch(
        fit_single_scdesign3_margin(bootstrap_counts, gene, pseudotime),
        error = function(error) NULL
      )
      if (is.null(bootstrap_fit)) {
        failed <- failed + 1L
        bootstrap_d[[b]] <- NA_real_
      } else {
        bootstrap_u <- randomized_nbi_pit(
          bootstrap_counts,
          bootstrap_fit$mu,
          bootstrap_fit$sigma
        )
        bootstrap_d[[b]] <- ks_stat_uniform(bootstrap_u)
      }
    }
    valid_d <- bootstrap_d[is.finite(bootstrap_d)]
    saveRDS(list(
      gene = gene, observed_u = observed_u, observed_d = observed_d,
      bootstrap_d = bootstrap_d, observed_seed = gene_seed,
      bootstrap_seed = gene_seed + 100000L, failed_fits = failed
    ), file.path(output_dir, "replicates", paste0(gene, ".rds")))
    if (failed > 0L || length(valid_d) != boot) {
      stop(sprintf("Failed bootstrap fits for %s: %d/%d; inspect saved replicates", gene, failed, boot))
    }
    exceedances <- sum(valid_d >= observed_d)
    mc_interval <- stats::binom.test(exceedances, boot)$conf.int
    row <- data.frame(
      gene = gene,
      cells = length(observed_counts),
      observed_ks_D = observed_d,
      naive_uniform_p = naive_p,
      refitted_bootstrap_p = (1 + exceedances) / (length(valid_d) + 1),
      exceedances = exceedances,
      tail_mc_lower = mc_interval[[1]],
      tail_mc_upper = mc_interval[[2]],
      bootstrap_q95 = as.numeric(stats::quantile(valid_d, 0.95, type = 1)),
      bootstrap_replicates = length(valid_d),
      failed_fits = failed,
      elapsed_seconds = proc.time()[["elapsed"]] - start,
      observed_seed = gene_seed,
      stringsAsFactors = FALSE
    )
    utils::write.csv(row, file.path(output_dir, "replicates", paste0(gene, ".csv")), row.names = FALSE)
    message(sprintf("Completed %s: B=%d, exceedances=%d", gene, boot, exceedances))
    row
  }
  rows <- parallel::mclapply(seq_along(genes), run_gene, mc.cores = workers, mc.set.seed = FALSE)
  if (any(vapply(rows, inherits, logical(1), "try-error"))) stop("A gene bootstrap failed")
  result <- do.call(rbind, rows)
  result$holm_refitted_p <- stats::p.adjust(result$refitted_bootstrap_p, method = "holm")
  result
}

if (Sys.getenv("SCDESIGN3_FUNCTIONS_ONLY", "0") != "1") {
data(example_sce, package = "scDesign3")
data(example_count, package = "scDesign3")
data(pseudotime, package = "scDesign3")

example_count_sce <- SingleCellExperiment(list(counts = as.matrix(example_count)))
colData(example_count_sce)$pseudotime <- pseudotime

diagnostics <- list(
  run_diagnostic("example_sce", example_sce, celltype_arg = "cell_type"),
  run_diagnostic("example_count", example_count_sce, celltype_arg = NULL)
)

metadata <- do.call(rbind, lapply(diagnostics, `[[`, "metadata"))
summary_df <- do.call(rbind, lapply(diagnostics, `[[`, "summary"))

bootstrap_b <- as.integer(Sys.getenv("SCDESIGN3_BOOTSTRAP_B", "99"))
bootstrap_workers <- as.integer(Sys.getenv("SCDESIGN3_WORKERS", "1"))
output_dir <- Sys.getenv("SCDESIGN3_RESULTS_DIR", "results/scdesign3_example")
bootstrap_genes <- head(rownames(example_count), 10)
refitted_bootstrap <- run_refitted_margin_bootstrap(
  counts = as.matrix(example_count),
  pseudotime = pseudotime,
  genes = bootstrap_genes,
  boot = bootstrap_b,
  workers = bootstrap_workers,
  output_dir = output_dir
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(
  refitted_bootstrap,
  file.path(output_dir, "pancreas_refitted_margin_bootstrap.csv"),
  row.names = FALSE
)

cat("\nscdesign3_example\n")
cat("dataset,cells,genes,cell_types,pseudotime_min,pseudotime_max\n")
for (i in seq_len(nrow(metadata))) {
  cat(sprintf(
    "%s,%d,%d,%s,%.3f,%.3f\n",
    metadata$dataset[i],
    metadata$cells[i],
    metadata$genes[i],
    metadata$cell_types[i],
    metadata$pseudotime_min[i],
    metadata$pseudotime_max[i]
  ))
}

cat("\ndataset,fit,median_gene_ks_D,max_gene_ks_D,min_naive_uniform_p,max_abs_pairwise_residual_corr,max_pseudotime_quartile_ks\n")
for (i in seq_len(nrow(summary_df))) {
  cat(sprintf(
    "%s,%s,%.3f,%.3f,%.3g,%.3f,%.3f\n",
    summary_df$dataset[i],
    summary_df$fit[i],
    summary_df$median_gene_ks_D[i],
    summary_df$max_gene_ks_D[i],
    summary_df$min_naive_uniform_p[i],
    summary_df$max_abs_pairwise_residual_corr[i],
    summary_df$max_pseudotime_quartile_ks[i]
  ))
}

cat("\nscdesign3_pancreas_refitted_bootstrap\n")
cat("gene,cells,observed_ks_D,naive_uniform_p,refitted_bootstrap_p,holm_refitted_p,bootstrap_q95,bootstrap_replicates,failed_fits,elapsed_seconds\n")
for (i in seq_len(nrow(refitted_bootstrap))) {
  cat(sprintf(
    "%s,%d,%.4f,%.4g,%.4g,%.4g,%.4f,%d,%d,%.3f\n",
    refitted_bootstrap$gene[i],
    refitted_bootstrap$cells[i],
    refitted_bootstrap$observed_ks_D[i],
    refitted_bootstrap$naive_uniform_p[i],
    refitted_bootstrap$refitted_bootstrap_p[i],
    refitted_bootstrap$holm_refitted_p[i],
    refitted_bootstrap$bootstrap_q95[i],
    refitted_bootstrap$bootstrap_replicates[i],
    refitted_bootstrap$failed_fits[i],
    refitted_bootstrap$elapsed_seconds[i]
  ))
}
}
