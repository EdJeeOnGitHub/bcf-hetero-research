# Heteroscedastic BCF Implementation Note

This fork adds an opt-in shared residual variance model for continuous-outcome
BCF:

```text
Y_i = mu_i + z_i tau_i + eps_i
eps_i ~ N(0, sigma_i^2)
sigma_i^2 = d(X_i)
```

The exported R entry point is `bcf_hetero()`. It calls the existing BCF sampler
with `variance_model = "shared"` and a `vartree` prior list. The returned object
contains posterior draw matrices `mu`, `tau`, `sigma2`, and `sigma`, where
`sigma = sqrt(sigma2)` on the original outcome scale. The scalar mean posterior
standard deviation used for legacy summaries is retained as `sigma_mean`.

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
against current squared residuals. When tree output is enabled, variance trees
are serialized alongside the control and moderator trees, and
`predict.bcf(..., type = "sigma2")` returns out-of-sample `sigma2` and `sigma`
draws.

Differences from `bayesm.HART`:

- scalar residual variance only;
- no full covariance or Cholesky off-diagonal trees;
- no DART/sparsity updates yet, though `vartree$sparse` is accepted and stored;
- variance prediction uses the same serialized-tree workflow as existing BCF
  prediction.
