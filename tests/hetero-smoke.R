library(bcf)

set.seed(1)
n <- 80
x <- matrix(rnorm(n * 3), n)
p <- pnorm(0.5 * x[, 1])
z <- rbinom(n, 1, p)
if (length(unique(z)) < 2) z[1:2] <- 0:1

mu <- x[, 1]
tau <- 0.5 + 0.25 * (x[, 2] > 0)
sigma2_true <- exp(0.4 * x[, 3])
y <- mu + z * tau + rnorm(n, sd = sqrt(sigma2_true))

tree_dir <- file.path(tempdir(), "bcf-hetero-smoke-trees")
dir.create(tree_dir, showWarnings = FALSE, recursive = TRUE)

fit <- bcf_hetero(
  y = y,
  z = z,
  x_control = x,
  x_moderate = x,
  pihat = p,
  x_variance = x,
  nburn = 3,
  nsim = 3,
  nthin = 1,
  n_chains = 1,
  n_threads = 1,
  ntree_control = 5,
  ntree_moderate = 3,
  vartree = list(num_trees = 3, nu = 5, lambda = 1, numcut = 20),
  no_output = FALSE,
  save_tree_directory = tree_dir,
  verbose = FALSE,
  use_muscale = FALSE,
  use_tauscale = FALSE
)

stopifnot(
  identical(dim(fit$mu), c(3L, as.integer(n))),
  identical(dim(fit$tau), c(3L, as.integer(n))),
  identical(dim(fit$sigma2), c(3L, as.integer(n))),
  identical(dim(fit$sigma), c(3L, as.integer(n))),
  isTRUE(all.equal(fit$vartree$lambda_tree, fit$vartree$lambda^(1 / fit$vartree$num_trees))),
  isTRUE(all.equal(fit$vartree$nu_tree, 2 / (1 - (1 - 2 / fit$vartree$nu)^(1 / fit$vartree$num_trees)))),
  all(is.finite(fit$sigma2)),
  all(fit$sigma2 > 0),
  isTRUE(all.equal(fit$sigma, sqrt(fit$sigma2), check.attributes = FALSE))
)

score <- bcf_hetero_cara_score(fit, alpha = 0.1)
stopifnot(
  length(score) == n,
  all(is.finite(score))
)

pred <- predict(
  fit,
  x_predict_control = x,
  x_predict_moderate = x,
  x_predict_variance = x,
  pi_pred = p,
  z_pred = z,
  save_tree_directory = tree_dir,
  type = "sigma2",
  n_cores = 1,
  verbose = FALSE
)

stopifnot(
  identical(dim(pred$sigma2), c(3L, as.integer(n))),
  identical(dim(pred$sigma), c(3L, as.integer(n))),
  all(is.finite(pred$sigma2)),
  all(pred$sigma2 > 0),
  isTRUE(all.equal(pred$sigma, sqrt(pred$sigma2), check.attributes = FALSE))
)

fit_ref <- bcf(
  y = y,
  z = z,
  x_control = x,
  x_moderate = x,
  pihat = p,
  nburn = 3,
  nsim = 3,
  nthin = 1,
  n_chains = 1,
  n_threads = 1,
  ntree_control = 5,
  ntree_moderate = 3,
  no_output = TRUE,
  verbose = FALSE,
  use_muscale = FALSE,
  use_tauscale = FALSE
)

diag <- bcf_hetero_diagnostics(
  fit,
  fit_reference = fit_ref,
  true_sigma2 = sigma2_true,
  alpha = 0.1,
  top_k = 10
)

stopifnot(
  length(diag$sigma2_mean) == n,
  is.finite(diag$sigma2_correlation) || is.na(diag$sigma2_correlation),
  length(diag$cara_score) == n,
  all(is.finite(diag$cara_score)),
  length(diag$rank_change) == n,
  diag$top_k == 10L,
  diag$top_k_overlap >= 0L,
  diag$top_k_overlap <= 10L,
  is.finite(diag$top_k_jaccard),
  diag$top_k_jaccard >= 0,
  diag$top_k_jaccard <= 1
)

cat("hetero-smoke OK\n")
