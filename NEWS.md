# bcf 2.0.2.9011 (research fork)

* Add optional joint shared/ratio variance scale proposals and paired variance
  leaf proposals, controlled by `joint_variance_every` and
  `paired_variance_every` (both disabled by default).
* Include optional variance-tree split changes and joint Gaussian mean-leaf
  refreshes from the intervening development builds.
* Include numerical posterior references for joint/paired variance moves,
  split-change reference checks, and saved-tree prediction replay checks.
* Regenerate Rcpp bindings to match the full sampler signature; this fixes
  stale bindings in the frozen source snapshot without changing sampler code.
* Document fork-specific installation, opt-in controls and validation limits.

## Experimental 2.0.2.9007

- Optional collapsed prognostic Cauchy amplitude removes redundant scale variables while preserving the marginal forest prior. Disabled by default; awaiting validation.

# bcf 2.0.2.9006 (experimental)

* Add optional `joint_mean_every` for an exact Gaussian refresh of all mean
  leaf coefficients at fixed partitions and scales. Default zero disables it.
* Record refresh counts and the largest coefficient block. Full posterior
  and saved-tree replay checks cover enabled and disabled settings.

# bcf 2.0.2.9005 (experimental)

* Add opt-in `global_mean_update`, a joint Gaussian refresh of the two
  ensemble intercept directions. It holds leaf contrasts fixed and preserves
  the current conditional prior and likelihood. Default FALSE.

# bcf 2.0.2.9004 (experimental)

* Add opt-in `paired_mean_update` for heteroskedastic fits. Joint Gaussian
  refreshes of prognostic and treatment tree leaves preserve the existing
  conditional prior and likelihood. The default remains disabled.
* Report successful and skipped paired updates. Analytical posterior and
  saved-tree replay checks exercise both settings before study validation.

# bcf 2.0.2.9003

* Correct the prognostic leaf-scale Gamma precision conditional to use its
  fixed base leaf variance. The old conditional included the previous precision
  in its residual statistic. Correct the corresponding optional non-half-Normal
  moderator branch as well.
* Add an almost-uninformative-likelihood check against the analytical log
  moment and tail probability implied by the documented half-Cauchy scale prior.

# bcf 2.0.2.9002

* Correct shared and treatment/control variance-ratio tree conditionals to
  include the observation precision `w` in squared residuals. Mean, scale,
  and variance updates now target the same documented weighted likelihood.
* Add weighted-variance regression checks against an analytical shared-variance
  posterior and a reference two-factor Gibbs sampler.

# bcf 2.0.2.9001

* Add `bcf_hetero()` for shared heteroscedastic residual variance using a scalar product-of-trees variance model.
* Store posterior draws of observation-level `sigma2` and add variance-aware CARA scoring/diagnostics helpers.
* Extend `predict.bcf()` to support variance prediction from saved heteroscedastic variance trees.
* Add a treatment-to-control variance-ratio model. Its log variance-ratio output
  is named `log_var_ratio`; the earlier development name `log_sigma_ratio` has
  been removed because it could be mistaken for a log standard-deviation ratio.
* Add fixed-scale fitting via `standardize = FALSE` and a joint prior-predictive
  SBC harness for prognostic, treatment-effect, homoscedastic variance, shared
  variance, and treatment/control variance-ratio recovery.

# bcf 2.0.2

### CRAN fixes

[Noah Greifer](https://github.com/ngreifer) updated the package source to reflect [two changes to the CRAN checks](https://www.tidyverse.org/blog/2023/03/cran-checks-compiled-code/)
that resulted in `bcf` being removed from CRAN in April 2023. Noah's updates:

1. Removed `sprintf()` from the C++ source code, as it is now deprecated, and 
2. Removed `CXX_STD = CXX11` from `src/Makevars` and `src/Makevars.win`, as C++11 is now a CRAN default.

### Serialization and performance updates

The prediction method introduced in the previous `bcf` version writes tree samples to text files, which can 
grow large if many samples are retained. Users concerned about the size of text file outputs 
may suppress writing to text files by specifying `no_output = TRUE` in the call to `bcf()`.

Sampling employs within-chain parallelism through `RcppParallel`, but `bcf` does not, 
for the time being, run multiple chains in parallel through R's high level `doParallel` interface.

## bcf 2.0.1

This implementation extends existing `bcf` functionality by:

- allowing for heteroskedastic errors
- automating multi-chain implementations
- providing a suite of convergence diagnostic functions via the `coda` package
- accelerating some underlying computations, resulting in shorter runtimes
- providing a function to predict treatment effects based on an existing model using new data

### Weights

The original version of `bcf` does not allow for weights, which are often used in practical applications to account for heteroskedasticity. Where the original BCF model was specified as:

y<sub>i</sub> &sim; N(&mu;(x<sub>i</sub>) + &tau;(x<sub>i</sub>) z<sub>i</sub>, &sigma;<sup>2</sup>),

which assumes that all outcomes y<sub>i</sub> have the same variance &sigma;<sup>2</sup>, in the extended version we can relax this assumption to allow for heteroskedasticity in y<sub>i</sub>:

y<sub>i</sub> &sim; N(&mu;(x<sub>i</sub>) + &tau;(x<sub>i</sub>) z<sub>i</sub>, &sigma;<sup>2</sup>/w<sub>i</sub>)

Incorporating weights impacts several parts of the code, including the computation of:

* sufficient statistics
* leaf node means
* leaf node means variance
* error variance (sigma)

### Automating multichain processing

In Bayesian analysis, it is useful to produce different runs of the same model -- with different starting values -- as a way of assessing convergence. If the different runs produce drastically different posterior distributions, it is a sign that the model has not converged fully.  In this version of `bcf` we have automated multichain processing and incorporated key MCMC diagnostics from the `coda` package, including effective sample sizes and the Gelman-Rubin statistic ("R hat").

### Within-chain parallelism

Finally, our implementation conducts some steps of the sampling procedure in parallel to maximize computational efficiency. Our testing shows that these enhancements have reduced runtimes by around 50%, across various experimental conditions.

### Implementing a prediction method

It is now possible to predict the treatment effect for a new set of units. Once users have produced a satisfactory `bcf` run (using training data), they can use this fitted `bcf` object to predict on a new set of test data. This is possible even with runs that have multiple chains.
