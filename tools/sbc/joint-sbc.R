#!/usr/bin/env Rscript

profile_config <- function(profile) {
  configs <- list(
    smoke = list(reps = 2L, n = 80L, nburn = 20L, nsim = 20L, nthin = 1L,
                 n_chains = 1L, ntree_control = 8L, ntree_moderate = 4L,
                 ntree_variance = 3L, workers = 1L),
    pilot = list(reps = 25L, n = 160L, nburn = 250L, nsim = 100L, nthin = 2L,
                 n_chains = 2L, ntree_control = 25L, ntree_moderate = 10L,
                 ntree_variance = 8L, workers = 4L),
    standard = list(reps = 250L, n = 240L, nburn = 500L, nsim = 200L, nthin = 2L,
                    n_chains = 2L, ntree_control = 40L, ntree_moderate = 15L,
                    ntree_variance = 16L, workers = 8L)
  )
  if (!profile %in% names(configs)) stop("profile must be smoke, pilot, or standard")
  c(configs[[profile]], list(
    profile = profile,
    sd_control = 1,
    sd_moderate = 0.5,
    base_control = 0.80,
    power_control = 2,
    base_moderate = 0.50,
    power_moderate = 3,
    sigma_nu = 5,
    sigma_lambda = 0.6,
    variance_nu = 5,
    variance_lambda = 0.6,
    variance_base = 0.80,
    variance_power = 2,
    variance_numcut = 50L,
    min_leaf = 5L
  ))
}

parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  out <- list(profile = "smoke", output = NULL, workers = NULL, reps = NULL,
              seed = 20260721L)
  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (grepl("^--[^=]+=", arg)) {
      key <- sub("^--([^=]+)=.*$", "\\1", arg)
      value <- sub("^--[^=]+=", "", arg)
    } else if (startsWith(arg, "--") && i < length(args)) {
      key <- sub("^--", "", arg)
      i <- i + 1L
      value <- args[[i]]
    } else {
      stop("unknown argument: ", arg)
    }
    if (!key %in% names(out)) stop("unknown option --", key)
    out[[key]] <- value
    i <- i + 1L
  }
  out$seed <- as.integer(out$seed)
  if (!is.null(out$workers)) out$workers <- as.integer(out$workers)
  if (!is.null(out$reps)) out$reps <- as.integer(out$reps)
  out
}

make_design <- function(n, seed) {
  set.seed(seed)
  x <- cbind(
    x1 = stats::runif(n, -1, 1),
    x2 = stats::rnorm(n),
    x3 = stats::runif(n, -1, 1)
  )
  pihat <- pmin(pmax(stats::pnorm(0.35 * x[, "x1"] - 0.20 * x[, "x2"]), 0.1), 0.9)
  z <- stats::rbinom(n, 1, pihat)
  if (sum(z) < 2L * 5L || sum(1L - z) < 2L * 5L) {
    ord <- order(pihat)
    z[ord[seq_len(10L)]] <- 0L
    z[rev(ord)[seq_len(10L)]] <- 1L
  }
  list(x = x, pihat = pihat, z = z)
}

make_prior_spec <- function(cfg, design) {
  x <- design$x
  cut_mean <- lapply(seq_len(ncol(x)), function(j) bcf:::.cp_quantile(x[, j]))
  cut_var <- lapply(seq_len(ncol(x)), function(j) {
    bcf:::.cp_quantile(x[, j], num = cfg$variance_numcut)
  })
  vartree <- bcf:::.parse_vartree(
    list(num_trees = cfg$ntree_variance, nu = cfg$variance_nu,
         lambda = cfg$variance_lambda, numcut = cfg$variance_numcut,
         base = cfg$variance_base, power = cfg$variance_power),
    yscale = rep(c(-1, 1), length.out = cfg$n), x_variance = x
  )
  list(cut_mean = cut_mean, cut_var = cut_var, vartree = vartree)
}

draw_truth <- function(model, cfg, design, spec, seed) {
  set.seed(seed)
  bcf:::bcfSbcPriorDraw(
    x_control = design$x,
    x_moderate = design$x,
    x_variance = design$x,
    z = design$z,
    cutpoints_control = spec$cut_mean,
    cutpoints_moderate = spec$cut_mean,
    cutpoints_variance = spec$cut_var,
    variance_model = model,
    ntree_control = cfg$ntree_control,
    ntree_moderate = cfg$ntree_moderate,
    ntree_variance = cfg$ntree_variance,
    sd_control = cfg$sd_control,
    sd_moderate = cfg$sd_moderate,
    base_control = cfg$base_control,
    power_control = cfg$power_control,
    base_moderate = cfg$base_moderate,
    power_moderate = cfg$power_moderate,
    variance_base = cfg$variance_base,
    variance_power = cfg$variance_power,
    sigma_nu = cfg$sigma_nu,
    sigma_lambda = cfg$sigma_lambda,
    variance_nu_tree = spec$vartree$nu_tree,
    variance_lambda_tree = spec$vartree$lambda_tree,
    min_leaf = cfg$min_leaf
  )
}

fit_model <- function(model, y, cfg, design, seed) {
  common <- list(
    y = y, z = design$z, x_control = design$x, x_moderate = design$x,
    pihat = design$pihat, x_variance = design$x,
    random_seed = seed, n_chains = cfg$n_chains, n_threads = 1L,
    nburn = cfg$nburn, nsim = cfg$nsim, nthin = cfg$nthin,
    ntree_control = cfg$ntree_control, ntree_moderate = cfg$ntree_moderate,
    sd_control = cfg$sd_control, sd_moderate = cfg$sd_moderate,
    base_control = cfg$base_control, power_control = cfg$power_control,
    base_moderate = cfg$base_moderate, power_moderate = cfg$power_moderate,
    nu = cfg$sigma_nu, lambda = cfg$sigma_lambda,
    include_pi = "none", use_muscale = FALSE, use_tauscale = FALSE,
    no_output = TRUE, verbose = FALSE, standardize = FALSE
  )
  suppressWarnings(if (identical(model, "homoscedastic")) {
    do.call(bcf::bcf, c(common, list(variance_model = model, vartree = NULL)))
  } else {
    do.call(bcf::bcf_hetero, c(common, list(
      variance_model = model,
      vartree = list(num_trees = cfg$ntree_variance, nu = cfg$variance_nu,
                     lambda = cfg$variance_lambda, numcut = cfg$variance_numcut,
                     base = cfg$variance_base, power = cfg$variance_power)
    )))
  })
}

randomized_rank <- function(draws, truth) {
  less <- sum(draws < truth)
  equal <- sum(draws == truth)
  less + if (equal > 0L) sample.int(equal + 1L, 1L) - 1L else 0L
}

chain_diagnostics <- function(draws, n_chains) {
  n <- length(draws)
  if (n_chains < 2L || n %% n_chains != 0L) {
    return(c(rhat = NA_real_, ess = unname(coda::effectiveSize(coda::mcmc(draws)))))
  }
  per_chain <- n %/% n_chains
  chains <- coda::mcmc.list(lapply(seq_len(n_chains), function(k) {
    ix <- ((k - 1L) * per_chain + 1L):(k * per_chain)
    coda::mcmc(draws[ix])
  }))
  rhat <- tryCatch(coda::gelman.diag(chains, autoburnin = FALSE,
                                    multivariate = FALSE)$psrf[1L, 1L],
                   error = function(e) NA_real_)
  c(rhat = unname(rhat), ess = unname(coda::effectiveSize(chains)))
}

rank_row <- function(replicate, scenario, estimand, draws, truth, cfg) {
  diag <- chain_diagnostics(draws, cfg$n_chains)
  qs <- stats::quantile(draws, c(0.05, 0.10, 0.25, 0.75, 0.90, 0.95), names = FALSE)
  data.frame(
    replicate = replicate, scenario = scenario, fit_model = scenario,
    estimand = estimand, truth = truth, posterior_mean = mean(draws),
    rank = randomized_rank(draws, truth), draws = length(draws),
    rhat = diag[["rhat"]], ess = diag[["ess"]],
    cover_50 = truth >= qs[3L] && truth <= qs[4L],
    cover_80 = truth >= qs[2L] && truth <= qs[5L],
    cover_90 = truth >= qs[1L] && truth <= qs[6L]
  )
}

rank_results <- function(fit, truth, scenario, replicate, point, cfg) {
  rows <- list(
    rank_row(replicate, scenario, "mean_mu", rowMeans(fit$mu), mean(truth$mu), cfg),
    rank_row(replicate, scenario, "point_mu", fit$mu[, point], truth$mu[point], cfg),
    rank_row(replicate, scenario, "mean_tau", rowMeans(fit$tau), mean(truth$tau), cfg),
    rank_row(replicate, scenario, "point_tau", fit$tau[, point], truth$tau[point], cfg)
  )
  if (identical(scenario, "homoscedastic")) {
    rows <- c(rows, list(rank_row(
      replicate, scenario, "sigma2", fit$sigma^2, truth$sigma0_2[1L], cfg
    )))
  } else {
    rows <- c(rows, list(
      rank_row(replicate, scenario, "mean_log_sigma0_2",
               rowMeans(log(fit$sigma0_2)), mean(log(truth$sigma0_2)), cfg),
      rank_row(replicate, scenario, "point_log_sigma0_2",
               log(fit$sigma0_2[, point]), log(truth$sigma0_2[point]), cfg)
    ))
    if (identical(scenario, "ratio")) {
      rows <- c(rows, list(
        rank_row(replicate, scenario, "mean_log_sigma1_2",
                 rowMeans(log(fit$sigma1_2)), mean(log(truth$sigma1_2)), cfg),
        rank_row(replicate, scenario, "point_log_sigma1_2",
                 log(fit$sigma1_2[, point]), log(truth$sigma1_2[point]), cfg),
        rank_row(replicate, scenario, "mean_log_var_ratio",
                 rowMeans(fit$log_var_ratio), mean(truth$log_var_ratio), cfg),
        rank_row(replicate, scenario, "point_log_var_ratio",
                 fit$log_var_ratio[, point], truth$log_var_ratio[point], cfg)
      ))
    }
  }
  do.call(rbind, rows)
}

coverage_rate <- function(draws, truth, prob = 0.9) {
  alpha <- (1 - prob) / 2
  ints <- apply(draws, 2L, stats::quantile,
                probs = c(alpha, 1 - alpha), names = FALSE)
  mean(truth >= ints[1L, ] & truth <= ints[2L, ])
}

log_mean_exp <- function(x) {
  m <- max(x)
  m + log(mean(exp(x - m)))
}

benchmark_row <- function(fit, truth, y_holdout, model, scenario, replicate, z) {
  mu_hat <- colMeans(fit$mu)
  tau_hat <- colMeans(fit$tau)
  if (identical(model, "homoscedastic")) {
    sigma0_hat <- sigma1_hat <- colMeans(fit$sigma2)
  } else if (identical(model, "shared")) {
    sigma0_hat <- sigma1_hat <- colMeans(fit$sigma2)
  } else {
    sigma0_hat <- colMeans(fit$sigma0_2)
    sigma1_hat <- colMeans(fit$sigma1_2)
  }
  sigma_draws <- if (identical(model, "ratio")) {
    fit$sigma0_2 * rep(1 - z, each = nrow(fit$sigma0_2)) +
      fit$sigma1_2 * rep(z, each = nrow(fit$sigma1_2))
  } else fit$sigma2
  mean_draws <- fit$mu + fit$tau * rep(z, each = nrow(fit$tau))
  log_score <- mean(vapply(seq_along(y_holdout), function(i) {
    log_mean_exp(stats::dnorm(y_holdout[i], mean_draws[, i],
                              sqrt(sigma_draws[, i]), log = TRUE))
  }, numeric(1)))
  data.frame(
    replicate = replicate, scenario = scenario, fit_model = model,
    mu_rmse = sqrt(mean((mu_hat - truth$mu)^2)),
    tau_rmse = sqrt(mean((tau_hat - truth$tau)^2)),
    log_variance_rmse = sqrt(mean(c(
      (log(sigma0_hat) - log(truth$sigma0_2))^2,
      (log(sigma1_hat) - log(truth$sigma1_2))^2
    ))),
    mu_coverage_90 = coverage_rate(fit$mu, truth$mu),
    tau_coverage_90 = coverage_rate(fit$tau, truth$tau),
    predictive_log_score = log_score
  )
}

run_replicate <- function(replicate, cfg, design, spec, seed) {
  scenarios <- c("homoscedastic", "shared", "ratio")
  rank_rows <- list()
  benchmark_rows <- list()
  for (s in seq_along(scenarios)) {
    scenario <- scenarios[[s]]
    rep_seed <- seed + replicate * 10000L + s * 1000L
    truth <- draw_truth(scenario, cfg, design, spec, rep_seed)
    set.seed(rep_seed + 1L)
    y <- truth$mean + stats::rnorm(cfg$n, sd = sqrt(truth$sigma2))
    y_holdout <- truth$mean + stats::rnorm(cfg$n, sd = sqrt(truth$sigma2))
    point <- sample.int(cfg$n, 1L)
    fit <- fit_model(scenario, y, cfg, design, rep_seed + 100L)
    rank_rows[[scenario]] <- rank_results(fit, truth, scenario, replicate, point, cfg)
    benchmark_rows[[paste0(scenario, "_correct")]] <- benchmark_row(
      fit, truth, y_holdout, scenario, scenario, replicate, design$z
    )
    if (!identical(scenario, "homoscedastic")) {
      base_fit <- fit_model("homoscedastic", y, cfg, design, rep_seed + 200L)
      benchmark_rows[[paste0(scenario, "_base")]] <- benchmark_row(
        base_fit, truth, y_holdout, "homoscedastic", scenario, replicate, design$z
      )
    }
  }
  list(ranks = do.call(rbind, rank_rows), benchmarks = do.call(rbind, benchmark_rows))
}

monte_carlo_interval <- function(n, probability, level = 0.99) {
  alpha <- (1 - level) / 2
  stats::qbinom(c(alpha, 1 - alpha), n, probability) / n
}

summarise_ranks <- function(ranks, profile) {
  ranks$included <- if (identical(profile, "smoke")) TRUE else
    !is.na(ranks$rhat) & ranks$rhat < 1.05 & ranks$ess >= 100
  groups <- split(ranks, interaction(ranks$scenario, ranks$estimand, drop = TRUE))
  out <- lapply(groups, function(df) {
    keep <- df[df$included, , drop = FALSE]
    n <- nrow(keep)
    u <- (keep$rank + 0.5) / (keep$draws + 1)
    sorted_u <- sort(u)
    ecdf_deviation <- if (n) {
      i <- seq_len(n)
      max(i / n - sorted_u, sorted_u - (i - 1) / n)
    } else NA_real_
    bins <- tabulate(pmin(10L, floor(u * 10) + 1L), nbins = 10L)
    p_chisq <- if (n >= 20L) suppressWarnings(stats::chisq.test(bins)$p.value) else NA_real_
    data.frame(
      scenario = df$scenario[1L], estimand = df$estimand[1L],
      replicates = nrow(df), included = n,
      excluded_rate = mean(!df$included),
      mean_rank_fraction = if (n) mean(u) else NA_real_,
      max_ecdf_deviation = ecdf_deviation,
      dkw_99_bound = if (n) sqrt(log(2 / 0.01) / (2 * n)) else NA_real_,
      chisq_p = p_chisq,
      coverage_50 = if (n) mean(keep$cover_50) else NA_real_,
      coverage_80 = if (n) mean(keep$cover_80) else NA_real_,
      coverage_90 = if (n) mean(keep$cover_90) else NA_real_,
      coverage_50_low = if (n) monte_carlo_interval(n, 0.5)[1L] else NA_real_,
      coverage_50_high = if (n) monte_carlo_interval(n, 0.5)[2L] else NA_real_,
      coverage_80_low = if (n) monte_carlo_interval(n, 0.8)[1L] else NA_real_,
      coverage_80_high = if (n) monte_carlo_interval(n, 0.8)[2L] else NA_real_,
      coverage_90_low = if (n) monte_carlo_interval(n, 0.9)[1L] else NA_real_,
      coverage_90_high = if (n) monte_carlo_interval(n, 0.9)[2L] else NA_real_
    )
  })
  summary <- do.call(rbind, out)
  summary$chisq_p_adjusted <- stats::p.adjust(summary$chisq_p, method = "BH")
  summary$calibrated <- summary$excluded_rate <= 0.05 &
    (is.na(summary$chisq_p_adjusted) | summary$chisq_p_adjusted >= 0.01) &
    summary$max_ecdf_deviation <= summary$dkw_99_bound &
    summary$coverage_50 >= summary$coverage_50_low & summary$coverage_50 <= summary$coverage_50_high &
    summary$coverage_80 >= summary$coverage_80_low & summary$coverage_80 <= summary$coverage_80_high &
    summary$coverage_90 >= summary$coverage_90_low & summary$coverage_90 <= summary$coverage_90_high
  summary$calibrated[is.na(summary$calibrated)] <- FALSE
  rownames(summary) <- NULL
  list(ranks = ranks, summary = summary)
}

summarise_benchmarks <- function(benchmarks) {
  metrics <- c("mu_rmse", "tau_rmse", "log_variance_rmse",
               "mu_coverage_90", "tau_coverage_90", "predictive_log_score")
  groups <- split(benchmarks, interaction(benchmarks$scenario,
                                           benchmarks$fit_model, drop = TRUE))
  do.call(rbind, lapply(groups, function(df) {
    values <- vapply(metrics, function(x) mean(df[[x]]), numeric(1))
    data.frame(scenario = df$scenario[1L], fit_model = df$fit_model[1L],
               as.list(values), row.names = NULL, check.names = FALSE)
  }))
}

pair_benchmarks <- function(benchmarks) {
  metrics <- c("mu_rmse", "tau_rmse", "log_variance_rmse",
               "mu_coverage_90", "tau_coverage_90", "predictive_log_score")
  base <- benchmarks[benchmarks$scenario != "homoscedastic" &
                       benchmarks$fit_model == "homoscedastic", , drop = FALSE]
  correct <- benchmarks[benchmarks$scenario != "homoscedastic" &
                          benchmarks$fit_model == benchmarks$scenario, , drop = FALSE]
  paired <- merge(correct, base, by = c("replicate", "scenario"),
                  suffixes = c("_hetero", "_base"), sort = TRUE)
  for (metric in metrics) {
    paired[[paste0(metric, "_difference")]] <-
      paired[[paste0(metric, "_hetero")]] - paired[[paste0(metric, "_base")]]
  }
  delta_names <- paste0(metrics, "_difference")
  groups <- split(paired, paired$scenario)
  summary <- do.call(rbind, lapply(groups, function(df) {
    means <- vapply(delta_names, function(x) mean(df[[x]]), numeric(1))
    ses <- vapply(delta_names, function(x) stats::sd(df[[x]]) / sqrt(nrow(df)), numeric(1))
    names(ses) <- paste0(delta_names, "_mcse")
    data.frame(scenario = df$scenario[1L], as.list(c(means, ses)),
               row.names = NULL, check.names = FALSE)
  }))
  list(replicates = paired, summary = summary)
}

plot_rank_histograms <- function(ranks, path) {
  keys <- unique(paste(ranks$scenario, ranks$estimand, sep = ": "))
  grDevices::png(path, width = 1800, height = 1200, res = 160)
  old <- graphics::par(mfrow = c(ceiling(length(keys) / 4), 4), mar = c(3, 3, 2, 1))
  on.exit({graphics::par(old); grDevices::dev.off()}, add = TRUE)
  for (key in keys) {
    ix <- paste(ranks$scenario, ranks$estimand, sep = ": ") == key & ranks$included
    u <- (ranks$rank[ix] + 0.5) / (ranks$draws[ix] + 1)
    graphics::hist(u, breaks = seq(0, 1, 0.1), main = key, xlab = "rank fraction",
                   col = "grey75", border = "white")
    if (length(u)) graphics::abline(h = length(u) / 10, col = "#B2182B", lwd = 2)
  }
}

plot_rank_ecdfs <- function(ranks, path) {
  keys <- unique(paste(ranks$scenario, ranks$estimand, sep = ": "))
  grDevices::png(path, width = 1800, height = 1200, res = 160)
  old <- graphics::par(mfrow = c(ceiling(length(keys) / 4), 4), mar = c(3, 3, 2, 1))
  on.exit({graphics::par(old); grDevices::dev.off()}, add = TRUE)
  for (key in keys) {
    ix <- paste(ranks$scenario, ranks$estimand, sep = ": ") == key & ranks$included
    u <- (ranks$rank[ix] + 0.5) / (ranks$draws[ix] + 1)
    if (length(u)) {
      graphics::plot(stats::ecdf(u), main = key, xlab = "rank fraction", ylab = "ECDF")
      graphics::abline(0, 1, col = "#2166AC", lwd = 2)
    } else graphics::plot.new()
  }
}

write_report <- function(path, cfg, rank_summary, benchmark_summary, elapsed) {
  failed <- rank_summary[!rank_summary$calibrated, c("scenario", "estimand"), drop = FALSE]
  lines <- c(
    "# Joint BCF simulation-based calibration",
    "",
    sprintf("- Profile: `%s`", cfg$profile),
    sprintf("- Replicates: %d", cfg$reps),
    sprintf("- Observations per replicate: %d", cfg$n),
    sprintf("- Elapsed minutes: %.1f", elapsed / 60),
    "- Correctly specified scenarios: homoscedastic, shared, ratio",
    "- Homoscedastic cross-fits on heteroscedastic draws are misspecification benchmarks, not SBC.",
    "",
    "## Calibration status",
    "",
    if (!nrow(failed)) "All monitored estimands passed." else
      paste("Failed:", paste(paste(failed$scenario, failed$estimand, sep = "/"), collapse = ", ")),
    "",
    paste0("See `rank-summary.csv`, `benchmark-summary.csv`, ",
           "`paired-benchmark-summary.csv`, `rank-histograms.png`, and `rank-ecdfs.png`.")
  )
  writeLines(lines, path)
}

run_joint_sbc <- function(cfg, output, seed) {
  dir.create(output, recursive = TRUE, showWarnings = FALSE)
  log_dir <- file.path(output, "logs")
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  design <- make_design(cfg$n, seed)
  spec <- make_prior_spec(cfg, design)
  started <- proc.time()[["elapsed"]]
  runner <- function(i) {
    con <- file(file.path(log_dir, sprintf("replicate-%04d.log", i)), open = "wt")
    sink(con, type = "output")
    sink(con, type = "message")
    on.exit({
      sink(type = "message")
      sink(type = "output")
      close(con)
    }, add = TRUE)
    run_replicate(i, cfg, design, spec, seed)
  }
  results <- if (cfg$workers > 1L && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(cfg$reps), runner, mc.cores = cfg$workers,
                       mc.preschedule = FALSE, mc.set.seed = FALSE)
  } else lapply(seq_len(cfg$reps), runner)
  failed_workers <- vapply(results, inherits, logical(1), what = "try-error")
  if (any(failed_workers)) {
    stop("SBC replicate workers failed: ", paste(which(failed_workers), collapse = ", "),
         "; see the corresponding log files")
  }
  elapsed <- proc.time()[["elapsed"]] - started
  ranks <- do.call(rbind, lapply(results, `[[`, "ranks"))
  benchmarks <- do.call(rbind, lapply(results, `[[`, "benchmarks"))
  calibrated <- summarise_ranks(ranks, cfg$profile)
  benchmark_summary <- summarise_benchmarks(benchmarks)
  paired_benchmarks <- pair_benchmarks(benchmarks)

  utils::write.csv(calibrated$ranks, file.path(output, "ranks.csv"), row.names = FALSE)
  utils::write.csv(calibrated$summary, file.path(output, "rank-summary.csv"), row.names = FALSE)
  utils::write.csv(benchmarks, file.path(output, "benchmarks.csv"), row.names = FALSE)
  utils::write.csv(benchmark_summary, file.path(output, "benchmark-summary.csv"), row.names = FALSE)
  utils::write.csv(paired_benchmarks$replicates,
                   file.path(output, "paired-benchmarks.csv"), row.names = FALSE)
  utils::write.csv(paired_benchmarks$summary,
                   file.path(output, "paired-benchmark-summary.csv"), row.names = FALSE)
  saveRDS(list(config = cfg, design = design, ranks = calibrated$ranks,
               rank_summary = calibrated$summary, benchmarks = benchmarks,
               benchmark_summary = benchmark_summary,
               paired_benchmarks = paired_benchmarks),
          file.path(output, "joint-sbc.rds"))
  plot_rank_histograms(calibrated$ranks, file.path(output, "rank-histograms.png"))
  plot_rank_ecdfs(calibrated$ranks, file.path(output, "rank-ecdfs.png"))
  write_report(file.path(output, "README.md"), cfg, calibrated$summary,
               benchmark_summary, elapsed)
  list(rank_summary = calibrated$summary, benchmark_summary = benchmark_summary,
       paired_benchmark_summary = paired_benchmarks$summary,
       elapsed = elapsed, output = output)
}

main <- function() {
  args <- parse_args()
  cfg <- profile_config(args$profile)
  if (!is.null(args$workers)) cfg$workers <- args$workers
  if (!is.null(args$reps)) cfg$reps <- args$reps
  output <- args$output
  if (is.null(output)) {
    output <- file.path("sbc-output", paste0(format(Sys.time(), "%Y%m%d-%H%M%S"),
                                              "-", cfg$profile))
  }
  cat(sprintf("joint-sbc profile=%s reps=%d workers=%d output=%s\n",
              cfg$profile, cfg$reps, cfg$workers, output))
  result <- run_joint_sbc(cfg, output, args$seed)
  print(result$rank_summary, row.names = FALSE)
  print(result$benchmark_summary, row.names = FALSE)
  print(result$paired_benchmark_summary, row.names = FALSE)
  if (identical(cfg$profile, "standard") && any(!result$rank_summary$calibrated)) {
    stop("one or more SBC calibration criteria failed; outputs were retained", call. = FALSE)
  }
  cat(sprintf("joint-sbc OK (%.1f minutes)\n", result$elapsed / 60))
}

if (sys.nframe() == 0L) main()
