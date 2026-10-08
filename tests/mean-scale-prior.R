library(bcf)

# A fixed prognostic stump has mu = scale * Normal(0, 1), where the documented
# marginal scale prior is half-Cauchy(1). Huge, tightly concentrated residual
# variance makes the likelihood effectively flat over this prior. This tests
# the implemented hyperparameter conditional against the model's analytical
# prior, rather than reproducing that conditional in another language.
n <- 40L
set.seed(7401)
x <- matrix(runif(n), ncol=1)
y <- rnorm(n)
samples <- numeric()
for(seed in c(7402,7403,7404)) {
  set.seed(seed)
  fit <- suppressWarnings(bcf(
    y=y, z=rep(0:1,n/2), x_control=x, x_moderate=x, pihat=rep(.5,n),
    nburn=1000, nsim=6000, n_chains=1, n_threads=1,
    ntree_control=1, ntree_moderate=1, sd_control=1, sd_moderate=1e-8,
    base_control=0, base_moderate=0, include_pi="none",
    collapsed_mu_scale=identical(Sys.getenv("BCF_COLLAPSED_MU_SCALE"),"true"),
    use_muscale=TRUE, use_tauscale=FALSE, standardize=FALSE,
    lambda=1e12, nu=1e6, no_output=TRUE, verbose=FALSE))
  samples <- c(samples,fit$mu[,1])
  rm(fit);gc()
}
expected_log_abs <- (digamma(.5)+log(2))/2
expected_tail <- integrate(function(theta)4/pi*pnorm(-1/tan(theta)),0,pi/2)$value
actual_log_abs <- mean(log(abs(samples)))
actual_tail <- mean(abs(samples)>1)
cat("Mean log absolute mu:",actual_log_abs,"expected:",expected_log_abs,"\n")
cat("Pr(abs(mu)>1):",actual_tail,"expected:",expected_tail,"\n")
stopifnot(all(is.finite(samples)),abs(actual_log_abs-expected_log_abs)<.12,
          abs(actual_tail-expected_tail)<.04)
cat("PASS: prognostic half-Cauchy scale prior matches analytical log moment and tail probability\n")
