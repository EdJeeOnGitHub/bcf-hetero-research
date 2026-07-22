library(bcf)

runner_candidates <- c(
  file.path("tools", "sbc", "joint-sbc.R"),
  file.path("..", "tools", "sbc", "joint-sbc.R"),
  file.path("..", "00_pkg_src", "bcf", "tools", "sbc", "joint-sbc.R")
)
runner <- runner_candidates[file.exists(runner_candidates)][1L]
if (is.na(runner)) stop("cannot locate tools/sbc/joint-sbc.R")
source(runner, local = TRUE)

## Fixed-scale inference must bypass every outcome transformation.
set.seed(20260721)
n <- 60L
x <- cbind(x1 = runif(n, -1, 1), x2 = rnorm(n))
pihat <- rep(0.5, n)
z <- rep(0:1, length.out = n)
y <- 3 + 2 * x[, 1] + rnorm(n, sd = 0.7)
fit_fixed <- suppressWarnings(bcf(
  y, z, x, x, pihat, nburn = 3, nsim = 3, n_chains = 1,
  n_threads = 1, ntree_control = 3, ntree_moderate = 2,
  sd_control = 1, sd_moderate = 0.5, nu = 5, lambda = 0.6,
  include_pi = "none", use_muscale = FALSE, use_tauscale = FALSE,
  no_output = TRUE, verbose = FALSE, standardize = FALSE
))
stopifnot(
  identical(fit_fixed$standardize, FALSE),
  identical(fit_fixed$muy, 0),
  identical(fit_fixed$sdy, 1)
)

## The prior generator is seeded, positive, and respects the variance-ratio API.
cfg <- profile_config("smoke")
design <- make_design(cfg$n, 20260721)
spec <- make_prior_spec(cfg, design)
set.seed(901)
draw_a <- draw_truth("ratio", cfg, design, spec, 901)
set.seed(901)
draw_b <- draw_truth("ratio", cfg, design, spec, 901)
stopifnot(
  identical(draw_a$mu, draw_b$mu),
  identical(draw_a$tau, draw_b$tau),
  identical(draw_a$sigma0_2, draw_b$sigma0_2),
  all(draw_a$sigma0_2 > 0),
  all(draw_a$sigma1_2 > 0),
  isTRUE(all.equal(
    draw_a$log_var_ratio,
    log(draw_a$sigma1_2 / draw_a$sigma0_2),
    check.attributes = FALSE
  )),
  all(draw_a$control$tree_sizes %% 2L == 1L),
  all(draw_a$variance$tree_sizes %% 2L == 1L)
)

## Basic prior-moment and topology checks guard the C++ generator itself.
moment_cfg <- cfg
moment_cfg$ntree_control <- moment_cfg$ntree_moderate <- moment_cfg$ntree_variance <- 1L
moment_cfg$base_control <- moment_cfg$base_moderate <- moment_cfg$variance_base <- 0
moment_spec <- make_prior_spec(moment_cfg, design)
moment_draws <- lapply(seq_len(300L), function(i) {
  draw_truth("shared", moment_cfg, design, moment_spec, 10000L + i)
})
control_leaves <- vapply(moment_draws, function(x) x$control$leaf_values[[1L]], numeric(1))
variance_leaves <- vapply(moment_draws, function(x) x$variance$leaf_values[[1L]], numeric(1))
expected_variance_mean <- moment_cfg$variance_nu * moment_cfg$variance_lambda /
  (moment_cfg$variance_nu - 2)
stopifnot(
  abs(mean(control_leaves)) < 0.15,
  abs(stats::sd(control_leaves) - moment_cfg$sd_control) < 0.15,
  abs(mean(variance_leaves) - expected_variance_mean) < 0.25
)

split_cfg <- moment_cfg
split_cfg$base_control <- 0.4
split_spec <- make_prior_spec(split_cfg, design)
root_split <- vapply(seq_len(400L), function(i) {
  draw_truth("homoscedastic", split_cfg, design, split_spec, 20000L + i)$control$tree_sizes[1L] > 1L
}, logical(1))
stopifnot(abs(mean(root_split) - split_cfg$base_control) < 0.1)

## Analytical inverse-chi-square rank calculation used by the SBC summaries.
set.seed(902)
nu <- 5
lambda <- 0.6
sigma2_true <- nu * lambda / rchisq(1, nu)
residuals <- rnorm(100, sd = sqrt(sigma2_true))
posterior_draws <- (nu * lambda + sum(residuals^2)) /
  rchisq(1000, nu + length(residuals))
rank <- randomized_rank(posterior_draws, sigma2_true)
stopifnot(rank >= 0L, rank <= length(posterior_draws))

## Exercise all three correct scenarios and both misspecified base cross-fits.
cfg$reps <- 2L
cfg$workers <- 1L
output <- file.path(tempdir(), "bcf-joint-sbc-smoke")
result <- run_joint_sbc(cfg, output, seed = 20260721)
stopifnot(
  all(c("homoscedastic", "shared", "ratio") %in% result$rank_summary$scenario),
  all(c("homoscedastic", "shared", "ratio") %in% result$benchmark_summary$fit_model),
  sum(result$benchmark_summary$fit_model == "homoscedastic") == 3L,
  all(file.exists(file.path(output, c(
    "joint-sbc.rds", "ranks.csv", "rank-summary.csv", "benchmarks.csv",
    "benchmark-summary.csv", "paired-benchmarks.csv",
    "paired-benchmark-summary.csv", "rank-histograms.png", "rank-ecdfs.png", "README.md"
  ))))
)

cat("joint-sbc-smoke OK\n")
