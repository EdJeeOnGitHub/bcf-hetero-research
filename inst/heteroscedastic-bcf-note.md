# Heteroscedastic BCF Implementation Note

This fork adds opt-in heteroscedastic residual variance models for
continuous-outcome BCF.

Shared mode:

```text
Y_i = mu_i + z_i tau_i + eps_i
eps_i ~ N(0, sigma_i^2)
sigma_i^2 = d(X_i)
```

Ratio mode:

```text
Y_i(0) = mu_i + eps0_i
Y_i(1) = mu_i + tau_i + eps1_i
eps0_i ~ N(0, sigma0_i^2)
eps1_i ~ N(0, sigma1_i^2)
sigma0_i^2 = d0(X_i)
sigma1_i^2 = sigma0_i^2 * r(X_i)
```

The exported R entry point is `bcf_hetero()`. It calls the existing BCF sampler
with `variance_model = "shared"` or `variance_model = "ratio"` and a `vartree`
prior list. The returned object contains posterior draw matrices `mu`, `tau`,
`sigma2`, and `sigma`, where `sigma = sqrt(sigma2)` on the original outcome
scale. In ratio mode it also contains `sigma0_2`, `sigma1_2`,
`log_var_ratio`, `sigma0`, and `sigma1`, where
`log_var_ratio = log(sigma1_2 / sigma0_2)`. The corresponding log standard-
deviation ratio is `log_var_ratio / 2`.

The variance function is a scalar product-of-trees model adapted from
Thomas Wiemann's `bayesm.HART` `varbart` machinery. Each variance tree has
positive inverse-chi-square leaves, and the residual variance prediction is the
product of tree predictions. Only the scalar diagonal path is implemented; the
modified-Cholesky `phitree` machinery is intentionally not included because BCF
needs a single residual variance per observation.

The user-facing `vartree$nu` and `vartree$lambda` are calibrated to per-tree
hyperparameters using the same product-of-trees transformation as the HART
reference:

```text
nu_tree = 2 / (1 - (1 - 2 / nu)^(1 / num_trees))
lambda_tree = lambda^(1 / num_trees)
```

At each MCMC iteration, the `mu` and `tau` tree updates use observation-specific
precision:

```text
omega_i = 1 / sigma_i^2
```

The existing BCF weighted sufficient-statistic path is used for these updates.
After the mean and treatment trees are updated, the variance trees are updated
against current squared residuals. In ratio mode, the baseline variance trees
are informed by both arms using the current treatment-arm ratio, and the ratio
trees are informed by treated residuals. When tree output is enabled, variance
trees are serialized alongside the control and moderator trees. Ratio mode also
serializes ratio trees. `predict.bcf(..., type = "sigma2")` returns
out-of-sample variance draws; ratio fits also return `sigma0_2`, `sigma1_2`,
and `log_var_ratio`.

Differences from `bayesm.HART`:

- scalar residual variance only;
- no full covariance or Cholesky off-diagonal trees;
- no DART/sparsity updates yet, though `vartree$sparse` is accepted and stored;
- variance prediction uses the same serialized-tree workflow as existing BCF
  prediction.

## Recovery Checks

`tests/hetero-smoke.R` is a fast API smoke test for shared and ratio modes. It
checks dimensions, finite positive variance draws, CARA scoring, diagnostics,
and saved-tree prediction.

`tests/hetero-recovery.R` is the sampler recovery harness. It repeatedly
simulates synthetic datasets with known mean, treatment-effect, and variance
functions, fits heteroscedastic BCF, and checks:

- correlation between posterior mean log-variance functions and truth;
- RMSE of recovered standard deviations;
- posterior rank of scalar truth summaries among posterior scalar draws;
- 90% posterior interval coverage for scalar variance summaries;
- extreme-rank rates for the SBC-style scalar diagnostics.

Run profiles:

```bash
Rscript tests/hetero-recovery.R --profile quick
Rscript tests/hetero-recovery.R --profile standard
BCF_RECOVERY_PROFILE=serious Rscript tests/hetero-recovery.R
```

The quick profile is intended for local iteration and gross failure detection.
The standard and serious profiles are the meaningful recovery checks. They are
SBC-style diagnostics rather than exact prior-SBC for the full nonparametric
BCF prior, because the data-generating functions are fixed synthetic truths
rather than draws from every BCF tree prior.

`tests/hetero-linear-posterior.R` compares BCF to an analytical conjugate
Bayesian linear-Gaussian posterior. The DGP is linear in `mu(x)` and `tau(x)`,
with known treatment-specific diagonal variances. The script computes the exact
normal posterior for the linear coefficients under a diffuse normal prior, then
compares BCF posterior draws of `mu_i` and `tau_i` to the analytical posterior
draws:

```bash
Rscript tests/hetero-linear-posterior.R --profile quick
Rscript tests/hetero-linear-posterior.R --profile standard
```

This is not an equality test between two identical Bayesian models. BCF has a
tree prior and estimates the variance functions, while the analytical linear
posterior conditions on known variances and uses a linear prior. The test is a
sanity check that, in a simple linear-Gaussian problem, BCF posterior estimands
land near the exact linear benchmark and have overlapping credible intervals.

## Joint Simulation-Based Calibration

`tools/sbc/joint-sbc.R` performs prior-predictive SBC of the prognostic,
treatment-effect, and variance components together. Its correctly specified
scenarios are homoscedastic, shared heteroscedastic, and treatment/control
variance ratio. The homoscedastic estimator is also fit to each heteroscedastic
draw as a misspecification benchmark; those cross-fit metrics are not treated
as SBC ranks.

The harness uses `standardize = FALSE`, fixed mean and treatment scales, fixed
variance hyperparameters, and an internal C++ generator that shares the
sampler's tree topology, cutpoint, minimum-leaf, normal-leaf, and scaled
inverse-chi-square laws. It writes replicate data, calibration and convergence
summaries, rank plots, benchmark metrics, and a Markdown report.

```bash
Rscript tools/sbc/joint-sbc.R --profile smoke
Rscript tools/sbc/joint-sbc.R --profile pilot
Rscript tools/sbc/joint-sbc.R --profile standard --workers 8
```
