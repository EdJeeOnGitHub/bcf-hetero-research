library(bcf)

parse_profile <- function(args = commandArgs(trailingOnly = TRUE)) {
  profile <- Sys.getenv("BCF_LINEAR_PROFILE", unset = "quick")
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
      n = 180L, nburn = 80L, nsim = 120L,
      ntree_control = 40L, ntree_moderate = 20L, vartree_num_trees = 12L,
      prior_sd = 10,
      max_mu_mean_abs_diff = 0.35,
      max_tau_mean_abs_diff = 0.35,
      min_mu_interval_overlap = 0.55,
      min_tau_interval_overlap = 0.55,
      min_bcf_truth_coverage = 0.70
    ),
    standard = list(
      n = 320L, nburn = 180L, nsim = 240L,
      ntree_control = 60L, ntree_moderate = 30L, vartree_num_trees = 20L,
      prior_sd = 10,
      max_mu_mean_abs_diff = 0.25,
      max_tau_mean_abs_diff = 0.25,
      min_mu_interval_overlap = 0.65,
      min_tau_interval_overlap = 0.65,
      min_bcf_truth_coverage = 0.80
    )
  )
  if (!profile %in% names(configs)) {
    stop("--profile must be one of: ", paste(names(configs), collapse = ", "),
         call. = FALSE)
  }
  cfg <- configs[[profile]]
  if (nzchar(Sys.getenv("BCF_LINEAR_N"))) {
    cfg$n <- as.integer(Sys.getenv("BCF_LINEAR_N"))
  }
  if (nzchar(Sys.getenv("BCF_LINEAR_NSIM"))) {
    cfg$nsim <- as.integer(Sys.getenv("BCF_LINEAR_NSIM"))
  }
  if (nzchar(Sys.getenv("BCF_LINEAR_NBURN"))) {
    cfg$nburn <- as.integer(Sys.getenv("BCF_LINEAR_NBURN"))
  }
  cfg
}

linear_design <- function(x, z) {
  cbind(
    intercept = 1,
    x1 = x[, "x1"],
    x2 = x[, "x2"],
    z = z,
    z_x1 = z * x[, "x1"],
    z_x2 = z * x[, "x2"]
  )
}

mu_design <- function(x) {
  cbind(
    intercept = 1,
    x1 = x[, "x1"],
    x2 = x[, "x2"],
    z = 0,
    z_x1 = 0,
    z_x2 = 0
  )
}

tau_design <- function(x) {
  cbind(
    intercept = 0,
    x1 = 0,
    x2 = 0,
    z = 1,
    z_x1 = x[, "x1"],
    z_x2 = x[, "x2"]
  )
}

posterior_linear_gaussian <- function(y, design, sigma2, prior_sd = 10) {
  stopifnot(length(y) == nrow(design), length(y) == length(sigma2))
  precision <- crossprod(design / sqrt(sigma2)) +
    diag(1 / prior_sd^2, ncol(design))
  cov <- solve(precision)
  mean <- cov %*% crossprod(design, y / sigma2)
  list(mean = drop(mean), cov = cov)
}

posterior_linear_draws <- function(post, ndraws, seed) {
  set.seed(seed)
  eig <- eigen(post$cov, symmetric = TRUE)
  transform <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = length(eig$values))
  z <- matrix(stats::rnorm(ndraws * length(post$mean)), nrow = ndraws)
  sweep(z %*% t(transform), 2, post$mean, "+")
}

draw_intervals <- function(draws, prob = 0.9) {
  alpha <- (1 - prob) / 2
  t(apply(draws, 2, stats::quantile, probs = c(alpha, 1 - alpha), names = FALSE))
}

interval_overlap_rate <- function(a, b) {
  lower <- pmax(a[, 1], b[, 1])
  upper <- pmin(a[, 2], b[, 2])
  mean(upper >= lower)
}

coverage_rate <- function(intervals, truth) {
  mean(truth >= intervals[, 1] & truth <= intervals[, 2])
}

make_data <- function(n, seed) {
  set.seed(seed)
  x <- cbind(
    x1 = stats::runif(n, -1, 1),
    x2 = stats::rnorm(n)
  )
  pihat <- stats::pnorm(0.20 * x[, "x1"] - 0.15 * x[, "x2"])
  z <- stats::rbinom(n, 1, pihat)
  if (length(unique(z)) < 2L) z[1:2] <- 0:1

  beta <- c(
    intercept = 0.15,
    x1 = 0.55,
    x2 = -0.30,
    z = 0.35,
    z_x1 = 0.25,
    z_x2 = -0.15
  )
  log_sigma0_2 <- -0.55 + 0.25 * x[, "x1"]
  log_ratio <- -0.30 + 0.20 * x[, "x2"]
  sigma0_2 <- exp(log_sigma0_2)
  sigma1_2 <- sigma0_2 * exp(log_ratio)
  sigma2_obs <- ifelse(z == 1L, sigma1_2, sigma0_2)
  y <- drop(linear_design(x, z) %*% beta) +
    stats::rnorm(n, sd = sqrt(sigma2_obs))

  list(
    x = x,
    z = z,
    pihat = pihat,
    y = y,
    beta = beta,
    sigma0_2 = sigma0_2,
    sigma1_2 = sigma1_2,
    sigma2_obs = sigma2_obs,
    mu_true = drop(mu_design(x) %*% beta),
    tau_true = drop(tau_design(x) %*% beta)
  )
}

fit_bcf_ratio <- function(dat, cfg) {
  bcf_hetero(
    y = dat$y,
    z = dat$z,
    x_control = dat$x,
    x_moderate = dat$x,
    x_variance = dat$x,
    pihat = dat$pihat,
    variance_model = "ratio",
    nburn = cfg$nburn,
    nsim = cfg$nsim,
    nthin = 1L,
    n_chains = 1L,
    n_threads = 1L,
    ntree_control = cfg$ntree_control,
    ntree_moderate = cfg$ntree_moderate,
    include_pi = "none",
    sd_control = 2,
    sd_moderate = 1,
    use_muscale = FALSE,
    use_tauscale = FALSE,
    vartree = list(
      num_trees = cfg$vartree_num_trees,
      nu = 5,
      lambda = 1,
      numcut = 50
    ),
    no_output = TRUE,
    verbose = FALSE
  )
}

main <- function() {
  profile <- parse_profile()
  cfg <- profile_config(profile)
  seed <- as.integer(Sys.getenv("BCF_LINEAR_SEED", unset = "20260703"))
  cat(sprintf(
    "hetero-linear-posterior profile=%s n=%d nburn=%d nsim=%d\n",
    profile, cfg$n, cfg$nburn, cfg$nsim
  ))

  dat <- make_data(cfg$n, seed)
  post <- posterior_linear_gaussian(
    y = dat$y,
    design = linear_design(dat$x, dat$z),
    sigma2 = dat$sigma2_obs,
    prior_sd = cfg$prior_sd
  )
  beta_draws <- posterior_linear_draws(post, cfg$nsim, seed + 1L)
  mu_lin_draws <- beta_draws %*% t(mu_design(dat$x))
  tau_lin_draws <- beta_draws %*% t(tau_design(dat$x))

  fit <- fit_bcf_ratio(dat, cfg)

  mu_lin_mean <- colMeans(mu_lin_draws)
  tau_lin_mean <- colMeans(tau_lin_draws)
  mu_bcf_mean <- colMeans(fit$mu)
  tau_bcf_mean <- colMeans(fit$tau)

  mu_lin_int <- draw_intervals(mu_lin_draws)
  tau_lin_int <- draw_intervals(tau_lin_draws)
  mu_bcf_int <- draw_intervals(fit$mu)
  tau_bcf_int <- draw_intervals(fit$tau)

  summary <- data.frame(
    profile = profile,
    n = cfg$n,
    mu_mean_abs_diff = mean(abs(mu_bcf_mean - mu_lin_mean)),
    tau_mean_abs_diff = mean(abs(tau_bcf_mean - tau_lin_mean)),
    mu_mean_cor = stats::cor(mu_bcf_mean, mu_lin_mean),
    tau_mean_cor = stats::cor(tau_bcf_mean, tau_lin_mean),
    mu_interval_overlap = interval_overlap_rate(mu_bcf_int, mu_lin_int),
    tau_interval_overlap = interval_overlap_rate(tau_bcf_int, tau_lin_int),
    mu_linear_truth_coverage = coverage_rate(mu_lin_int, dat$mu_true),
    tau_linear_truth_coverage = coverage_rate(tau_lin_int, dat$tau_true),
    mu_bcf_truth_coverage = coverage_rate(mu_bcf_int, dat$mu_true),
    tau_bcf_truth_coverage = coverage_rate(tau_bcf_int, dat$tau_true),
    sigma0_cor = stats::cor(log(colMeans(fit$sigma0_2)), log(dat$sigma0_2)),
    sigma1_cor = stats::cor(log(colMeans(fit$sigma1_2)), log(dat$sigma1_2))
  )
  print(summary, row.names = FALSE)

  checks <- c(
    summary$mu_mean_abs_diff <= cfg$max_mu_mean_abs_diff,
    summary$tau_mean_abs_diff <= cfg$max_tau_mean_abs_diff,
    summary$mu_interval_overlap >= cfg$min_mu_interval_overlap,
    summary$tau_interval_overlap >= cfg$min_tau_interval_overlap,
    summary$mu_bcf_truth_coverage >= cfg$min_bcf_truth_coverage,
    summary$tau_bcf_truth_coverage >= cfg$min_bcf_truth_coverage
  )
  names(checks) <- c(
    "mu posterior mean distance",
    "tau posterior mean distance",
    "mu interval overlap",
    "tau interval overlap",
    "BCF mu truth coverage",
    "BCF tau truth coverage"
  )

  failed <- names(checks)[!checks]
  if (length(failed) > 0L) {
    stop(
      "Linear analytical posterior comparison failed for profile '", profile, "': ",
      paste(failed, collapse = ", "),
      call. = FALSE
    )
  }
  cat("hetero-linear-posterior OK\n")
}

main()
