library(bcf)

parse_profile <- function(args = commandArgs(trailingOnly = TRUE)) {
  profile <- Sys.getenv("BCF_RECOVERY_PROFILE", unset = "quick")
  for (i in seq_along(args)) {
    if (identical(args[[i]], "--profile") && i < length(args)) {
      profile <- args[[i + 1L]]
    } else if (startsWith(args[[i]], "--profile=")) {
      profile <- sub("^--profile=", "", args[[i]])
    }
  }
  profile
}

profile_config <- function(profile) {
  configs <- list(
    quick = list(
      reps = 2L, n = 140L, nburn = 30L, nsim = 40L,
      ntree_control = 20L, ntree_moderate = 8L, vartree_num_trees = 8L,
      min_shared_log_sigma_cor = 0.25,
      min_ratio_log_sigma0_cor = 0.15,
      min_ratio_log_sigma1_cor = 0.15,
      min_scalar_coverage = 0.25,
      max_extreme_rank_rate = 0.75
    ),
    standard = list(
      reps = 8L, n = 240L, nburn = 120L, nsim = 160L,
      ntree_control = 40L, ntree_moderate = 15L, vartree_num_trees = 16L,
      min_shared_log_sigma_cor = 0.45,
      min_ratio_log_sigma0_cor = 0.30,
      min_ratio_log_sigma1_cor = 0.30,
      min_scalar_coverage = 0.50,
      max_extreme_rank_rate = 0.50
    ),
    serious = list(
      reps = 24L, n = 360L, nburn = 250L, nsim = 300L,
      ntree_control = 60L, ntree_moderate = 20L, vartree_num_trees = 24L,
      min_shared_log_sigma_cor = 0.55,
      min_ratio_log_sigma0_cor = 0.40,
      min_ratio_log_sigma1_cor = 0.40,
      min_scalar_coverage = 0.65,
      max_extreme_rank_rate = 0.35
    )
  )
  if (!profile %in% names(configs)) {
    stop("--profile must be one of: ", paste(names(configs), collapse = ", "),
         call. = FALSE)
  }
  cfg <- configs[[profile]]
  if (nzchar(Sys.getenv("BCF_RECOVERY_REPS"))) {
    cfg$reps <- as.integer(Sys.getenv("BCF_RECOVERY_REPS"))
  }
  if (nzchar(Sys.getenv("BCF_RECOVERY_N"))) {
    cfg$n <- as.integer(Sys.getenv("BCF_RECOVERY_N"))
  }
  if (nzchar(Sys.getenv("BCF_RECOVERY_NSIM"))) {
    cfg$nsim <- as.integer(Sys.getenv("BCF_RECOVERY_NSIM"))
  }
  if (nzchar(Sys.getenv("BCF_RECOVERY_NBURN"))) {
    cfg$nburn <- as.integer(Sys.getenv("BCF_RECOVERY_NBURN"))
  }
  cfg
}

cor_complete <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3L || stats::sd(x[ok]) == 0 || stats::sd(y[ok]) == 0) {
    return(NA_real_)
  }
  stats::cor(x[ok], y[ok])
}

rmse <- function(x, y) sqrt(mean((x - y)^2))

draw_rank <- function(draws, truth) {
  mean(draws < truth)
}

central_interval_covers <- function(draws, truth, prob = 0.9) {
  alpha <- (1 - prob) / 2
  q <- stats::quantile(draws, probs = c(alpha, 1 - alpha), names = FALSE)
  truth >= q[[1L]] && truth <= q[[2L]]
}

make_design <- function(n) {
  cbind(
    x1 = stats::runif(n, -1, 1),
    x2 = stats::rnorm(n),
    x3 = stats::runif(n, -1, 1)
  )
}

make_assignment <- function(pihat) {
  z <- stats::rbinom(length(pihat), 1, pihat)
  if (length(unique(z)) < 2L) {
    z[seq_len(min(2L, length(pihat)))] <- 0:1
  }
  z
}

fit_args <- function(cfg) {
  list(
    nburn = cfg$nburn,
    nsim = cfg$nsim,
    nthin = 1L,
    n_chains = 1L,
    n_threads = 1L,
    ntree_control = cfg$ntree_control,
    ntree_moderate = cfg$ntree_moderate,
    vartree = list(
      num_trees = cfg$vartree_num_trees,
      nu = 5,
      lambda = 1,
      numcut = 50
    ),
    no_output = TRUE,
    verbose = FALSE,
    use_muscale = FALSE,
    use_tauscale = FALSE
  )
}

fit_hetero <- function(y, z, x, pihat, variance_model, cfg) {
  do.call(
    bcf_hetero,
    c(
      list(
        y = y,
        z = z,
        x_control = x,
        x_moderate = x,
        x_variance = x,
        pihat = pihat,
        variance_model = variance_model
      ),
      fit_args(cfg)
    )
  )
}

run_shared_rep <- function(seed, cfg) {
  set.seed(seed)
  x <- make_design(cfg$n)
  pihat <- stats::pnorm(0.35 * x[, "x1"] - 0.20 * x[, "x2"])
  z <- make_assignment(pihat)
  mu <- 0.45 * x[, "x1"] - 0.25 * x[, "x2"]
  tau <- 0.35 + 0.25 * (x[, "x1"] > 0)
  log_sigma2 <- -0.35 + 0.95 * (x[, "x3"] > 0) + 0.15 * x[, "x2"]
  sigma2 <- exp(log_sigma2)
  y <- mu + z * tau + stats::rnorm(cfg$n, sd = sqrt(sigma2))

  fit <- fit_hetero(y, z, x, pihat, "shared", cfg)
  log_sigma2_draws <- log(pmax(fit$sigma2, .Machine$double.eps))
  mean_log_sigma2_draws <- rowMeans(log_sigma2_draws)
  tau_draws <- fit$tau

  data.frame(
    scenario = "shared",
    seed = seed,
    log_sigma_cor = cor_complete(colMeans(log_sigma2_draws), log_sigma2),
    log_sigma0_cor = NA_real_,
    log_sigma1_cor = NA_real_,
    log_ratio_cor = NA_real_,
    sigma_rmse = rmse(sqrt(colMeans(fit$sigma2)), sqrt(sigma2)),
    sigma0_rmse = NA_real_,
    sigma1_rmse = NA_real_,
    tau_cor = cor_complete(colMeans(tau_draws), tau),
    scalar = "mean_log_sigma2",
    scalar_truth = mean(log_sigma2),
    scalar_post_mean = mean(mean_log_sigma2_draws),
    scalar_rank = draw_rank(mean_log_sigma2_draws, mean(log_sigma2)),
    scalar_covered_90 = central_interval_covers(mean_log_sigma2_draws, mean(log_sigma2))
  )
}

run_ratio_rep <- function(seed, cfg) {
  set.seed(seed)
  x <- make_design(cfg$n)
  pihat <- stats::pnorm(0.35 * x[, "x1"] - 0.20 * x[, "x2"])
  z <- make_assignment(pihat)
  mu <- 0.40 * x[, "x1"] - 0.30 * x[, "x2"]
  tau <- 0.30 + 0.20 * (x[, "x1"] > 0)
  log_sigma0_2 <- -0.45 + 0.90 * (x[, "x3"] > 0) + 0.10 * x[, "x2"]
  log_ratio <- -0.50 + 0.75 * (x[, "x2"] > 0)
  log_sigma1_2 <- log_sigma0_2 + log_ratio
  sigma0_2 <- exp(log_sigma0_2)
  sigma1_2 <- exp(log_sigma1_2)
  y <- mu + z * tau + stats::rnorm(
    cfg$n,
    sd = sqrt(ifelse(z == 1L, sigma1_2, sigma0_2))
  )

  fit <- fit_hetero(y, z, x, pihat, "ratio", cfg)
  log_sigma0_2_draws <- log(pmax(fit$sigma0_2, .Machine$double.eps))
  log_sigma1_2_draws <- log(pmax(fit$sigma1_2, .Machine$double.eps))
  log_ratio_draws <- fit$log_var_ratio

  data.frame(
    scenario = "ratio",
    seed = seed,
    log_sigma_cor = NA_real_,
    log_sigma0_cor = cor_complete(colMeans(log_sigma0_2_draws), log_sigma0_2),
    log_sigma1_cor = cor_complete(colMeans(log_sigma1_2_draws), log_sigma1_2),
    log_ratio_cor = cor_complete(colMeans(log_ratio_draws), log_ratio),
    sigma_rmse = NA_real_,
    sigma0_rmse = rmse(sqrt(colMeans(fit$sigma0_2)), sqrt(sigma0_2)),
    sigma1_rmse = rmse(sqrt(colMeans(fit$sigma1_2)), sqrt(sigma1_2)),
    tau_cor = cor_complete(colMeans(fit$tau), tau),
    scalar = "mean_log_ratio",
    scalar_truth = mean(log_ratio),
    scalar_post_mean = mean(rowMeans(log_ratio_draws)),
    scalar_rank = draw_rank(rowMeans(log_ratio_draws), mean(log_ratio)),
    scalar_covered_90 = central_interval_covers(rowMeans(log_ratio_draws), mean(log_ratio))
  )
}

summarise_recovery <- function(results) {
  split_results <- split(results, results$scenario)
  do.call(rbind, lapply(split_results, function(df) {
    numeric_cols <- names(df)[vapply(df, is.numeric, logical(1))]
    numeric_cols <- setdiff(numeric_cols, "seed")
    out <- as.data.frame(as.list(vapply(df[numeric_cols], mean, numeric(1), na.rm = TRUE)))
    out$scenario <- df$scenario[[1L]]
    out$reps <- nrow(df)
    out$scalar_coverage_90 <- mean(df$scalar_covered_90, na.rm = TRUE)
    out$scalar_extreme_rank_rate <- mean(df$scalar_rank <= 0.05 | df$scalar_rank >= 0.95,
                                         na.rm = TRUE)
    out[, c("scenario", "reps", setdiff(names(out), c("scenario", "reps")))]
  }))
}

assert_recovery <- function(summary, cfg, profile) {
  shared <- summary[summary$scenario == "shared", , drop = FALSE]
  ratio <- summary[summary$scenario == "ratio", , drop = FALSE]
  if (nrow(shared) != 1L || nrow(ratio) != 1L) {
    stop("Expected one shared and one ratio summary row.", call. = FALSE)
  }

  checks <- c(
    shared$log_sigma_cor >= cfg$min_shared_log_sigma_cor,
    ratio$log_sigma0_cor >= cfg$min_ratio_log_sigma0_cor,
    ratio$log_sigma1_cor >= cfg$min_ratio_log_sigma1_cor,
    shared$scalar_coverage_90 >= cfg$min_scalar_coverage,
    ratio$scalar_coverage_90 >= cfg$min_scalar_coverage,
    shared$scalar_extreme_rank_rate <= cfg$max_extreme_rank_rate,
    ratio$scalar_extreme_rank_rate <= cfg$max_extreme_rank_rate
  )
  names(checks) <- c(
    "shared log-sigma correlation",
    "ratio log-sigma0 correlation",
    "ratio log-sigma1 correlation",
    "shared scalar 90% coverage",
    "ratio scalar 90% coverage",
    "shared scalar extreme rank rate",
    "ratio scalar extreme rank rate"
  )

  failed <- names(checks)[!checks]
  if (length(failed) > 0L) {
    print(summary, row.names = FALSE)
    stop(
      "Heteroscedastic BCF recovery checks failed for profile '", profile, "': ",
      paste(failed, collapse = ", "),
      call. = FALSE
    )
  }
}

main <- function() {
  profile <- parse_profile()
  cfg <- profile_config(profile)
  cat(sprintf(
    "hetero-recovery profile=%s reps=%d n=%d nburn=%d nsim=%d\n",
    profile, cfg$reps, cfg$n, cfg$nburn, cfg$nsim
  ))

  base_seed <- as.integer(Sys.getenv("BCF_RECOVERY_SEED", unset = "20260702"))
  shared <- do.call(rbind, lapply(seq_len(cfg$reps), function(i) {
    run_shared_rep(base_seed + i, cfg)
  }))
  ratio <- do.call(rbind, lapply(seq_len(cfg$reps), function(i) {
    run_ratio_rep(base_seed + 10000L + i, cfg)
  }))

  results <- rbind(shared, ratio)
  summary <- summarise_recovery(results)
  print(summary, row.names = FALSE)
  assert_recovery(summary, cfg, profile)
  cat("hetero-recovery OK\n")
}

main()
