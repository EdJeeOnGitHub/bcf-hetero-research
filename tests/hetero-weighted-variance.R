library(bcf)

# With mean functions pinned near zero and unsplittable variance stumps,
# the shared variance posterior is analytically inverse chi-square. This
# exercises the C++ updates under non-unit likelihood precisions, including
# weights below one, rather than comparing two implementations of a formula.
set.seed(7301)
n <- 180L
x <- matrix(runif(n), ncol=1)
z <- rep(0:1, length.out=n)
w <- rep(c(0.25, 1, 9), length.out=n)
y <- rnorm(n, sd=0.8)
nu <- 8
lambda <- 0.6
run_fit <- function(model, seed) suppressWarnings(bcf_hetero(
  y=y, z=z, x_control=x, x_moderate=x, x_variance=x,
  pihat=rep(0.5,n), w=w, random_seed=seed, n_chains=1, n_threads=1,
  nburn=1000, nsim=6000, ntree_control=1, ntree_moderate=1,
  sd_control=1e-8, sd_moderate=1e-8,
  base_control=0, base_moderate=0, include_pi="none",
  use_muscale=FALSE, use_tauscale=FALSE, standardize=FALSE,
  vartree=list(num_trees=1, nu=nu, lambda=lambda, base=0),
  variance_model=model, no_output=TRUE, verbose=FALSE
))

shared <- run_fit("shared", 7302)
rss <- sum(w*y^2)
expected <- (nu*lambda+rss)/(nu+n-2)
expected_quantiles <- (nu*lambda+rss)/qchisq(c(.9,.5,.1),nu+n)
stopifnot(
  max(abs(shared$mu))<1e-5,
  abs(mean(shared$sigma0_2[,1])/expected-1)<.05,
  max(abs(unname(quantile(shared$sigma0_2[,1],c(.1,.5,.9)))/expected_quantiles-1))<.05
)

# The ratio model's two scalar variance factors have a simple independent
# reference Gibbs sampler. Compare both potential-state variance marginals;
# both the baseline-tree and ratio-tree weighted residual paths matter here.
set.seed(7303)
rss0 <- sum(w[z==0]*y[z==0]^2)
rss1 <- sum(w[z==1]*y[z==1]^2)
v0 <- ratio <- 1
ref <- matrix(NA_real_, 60000, 2)
for(i in seq_len(nrow(ref))) {
  v0 <- (nu*lambda + rss0 + rss1/ratio)/rchisq(1,nu+n)
  ratio <- (nu*lambda + rss1/v0)/rchisq(1,nu+sum(z==1))
  ref[i,] <- c(v0,v0*ratio)
}
ref <- ref[-seq_len(10000),,drop=FALSE]
fit_ratio <- run_fit("ratio",7304)
actual <- cbind(fit_ratio$sigma0_2[,1],fit_ratio$sigma1_2[,1])
stopifnot(
  max(abs(colMeans(actual)/colMeans(ref)-1))<.08,
  max(abs(apply(actual,2,quantile,c(.1,.5,.9))/apply(ref,2,quantile,c(.1,.5,.9))-1))<.08
)
cat("PASS: weighted shared and ratio variance posteriors match analytical/reference distributions\n")
