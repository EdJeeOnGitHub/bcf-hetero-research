# Heteroscedastic Bayesian Causal Forests

This research fork of [jaredsmurray/bcf](https://github.com/jaredsmurray/bcf)
adds covariate-dependent residual variance, treatment/control variance ratios,
and optional MCMC moves for continuous-outcome Bayesian causal forests.
The package name remains `bcf`. The prepared release is **2.0.2.9011**.
It is a development fork, not the CRAN release.

## Install this fork

Install the tagged source once the release has been pushed:

```r
install.packages("remotes")
remotes::install_github("EdJeeOnGitHub/bcf-hetero-private@v2.0.2.9011")
```

The repository currently requires access while private. Changing its visibility
is a separate owner action; installation will become anonymous once public.
The URL above must be updated if the repository is renamed.

For a local checkout:

```sh
R CMD build --no-build-vignettes --no-manual .
R CMD INSTALL bcf_2.0.2.9011.tar.gz
```

Source installation needs a C++ toolchain and the dependencies declared in
`DESCRIPTION`; Windows users need Rtools. Installing CRAN `bcf` does not install
these fork extensions.

## Model and sampler options

Use `variance_model = "shared"` for covariate-dependent residual variance, or
`variance_model = "ratio"` for separate treatment/control variances linked by a
ratio. Ratio fits return `sigma0_2`, `sigma1_2`, and
`log_var_ratio = log(sigma1_2 / sigma0_2)`. See
[the implementation note](inst/heteroscedastic-bcf-note.md) for the likelihood,
weights, priors, prediction support and reference implementation attribution.
Shared fits return `sigma2`; `sigma2`, `sigma0_2`, and `sigma1_2` are in squared
outcome units. The log variance ratio is dimensionless, and its half is the
log standard-deviation ratio. With observation precision weights `w`, the
likelihood variance is the corresponding variance output divided by `w`.

The following optional moves leave the specified target posterior unchanged:

| Argument | Purpose | Package default |
|---|---|---|
| `joint_mean_every` | Joint Gaussian refresh of mean leaf coefficients | `0` (off) |
| `variance_split_change` | Change existing variance-tree split rules | `FALSE` |
| `joint_variance_every` | Joint shared/ratio variance scale proposals | `0` (off) |
| `paired_variance_every` | Paired shared/ratio variance leaf proposals | `0` (off) |

For ratio fits, a common configuration is:

```r
# Add these arguments to bcf(..., variance_model = "ratio", ...):
joint_mean_every = 5L
variance_split_change = TRUE
joint_variance_every = 1L
paired_variance_every = 1L
```

Joint and paired variance moves require the ratio model. These controls change
sampling, not the prior. Choose prior scales, chain lengths and convergence
criteria for the analysis; this configuration is not a convergence guarantee.
Other experimental controls remain opt-in; their presence is not an endorsement
for every application.

## Validation

The release preparation checks include analytical mean-conditionals for fixed
and half-Cauchy scales, numerical joint/paired variance posterior references,
variance split proposal references, and saved-tree prediction replay. See
[release validation](docs/release-9011.md) for commands and results.

[Historical simulation-based calibration](inst/sbc-results/standard-20260721/README.md)
includes reported failures as well as passes. It predates this release and must
not be interpreted as a fresh calibration of the new sampling configuration.
Implementation checks and synthetic recovery do not establish convergence or
inferential validity for any particular real dataset.

## Attribution and license

This fork preserves upstream BCF author attribution in `DESCRIPTION` and
`inst/CITATION`. The original method is Hahn, Murray and Carvalho (2020).
The scalar product-of-trees variance implementation is adapted from Thomas
Wiemann's `bayesm.HART` machinery, as described in the implementation note.
The package is licensed under **GPL-3**, as declared in `DESCRIPTION`.

## Upstream introduction

# Bayesian Causal Forests

Welcome to the `bcf` site! This page provides hands-on examples of how to conduct Bayesian causal forest (BCF) analyses. You can find methodological details on the underlying modeling approach in the original BCF paper: [Hahn, Murray, and Carvalho 2020](https://projecteuclid.org/journals/bayesian-analysis/volume-15/issue-3/Bayesian-Regression-Tree-Models-for-Causal-Inference--Regularization-Confounding/10.1214/19-BA1195.full).

## Why BCF?

BCF is a cutting-edge model for causal inference that builds on Bayesian Additive Regression Trees (BART, [Chipman, George, and McCulloch 2010](https://projecteuclid.org/euclid.aoas/1273584455)). BART and BCF both combine Bayesian regularization with regression trees to provide a highly flexible response surface that, thanks to regularization from prior distributions, does not overfit to the training data. BCF extends BART's flexibility by specifying different models for relationships between (1) covariates and the outcome and (2) covariates and the treatment effect, and regularizing the treatment effect directly.

BCF performs remarkably well in simulation and has led the pack at recent rigorous causal inference competitions, such as those held at the Atlantic Causal Inference Conference (see, for example, [Dorie et al. 2019](https://projecteuclid.org/euclid.ss/1555056030)).

## Getting Started

If you are just getting started with `bcf`, we recommend beginning with the tutorial vignettes.

## Heteroscedastic residual variance

`bcf_hetero()` supports a shared variance function and a treatment-to-control
variance-ratio model. Ratio fits and predictions return `sigma0_2`,
`sigma1_2`, and `log_var_ratio`, with
`log_var_ratio = log(sigma1_2 / sigma0_2)`. Thus the log standard-deviation
ratio is `log_var_ratio / 2`.

## Simulation-based calibration

The joint SBC harness draws prognostic, treatment-effect, and variance trees
from their exact fixed-scale priors, simulates outcomes, and refits the
homoscedastic, shared-variance, and variance-ratio models. It also compares the
homoscedastic estimator on the same heteroscedastic draws as a deliberately
misspecified benchmark.

```sh
Rscript tools/sbc/joint-sbc.R --profile smoke
Rscript tools/sbc/joint-sbc.R --profile pilot
Rscript tools/sbc/joint-sbc.R --profile standard --workers 8
```

The standard profile runs 250 replicates and writes rank data, coverage and
convergence summaries, benchmark metrics, plots, and a Markdown report beneath
`sbc-output/`. Exact prior-predictive calibration uses `standardize = FALSE`
and explicitly fixed prior scales; ordinary fits continue to standardize by
default.

## Installation

This package requires compilation, so make sure you have Rtools properly installed if you are on Windows -- see [this site](https://cran.r-project.org/bin/windows/Rtools/) for details.

Install the latest release from CRAN:

```{r}
install.packages("bcf")
```

Install the upstream development version (without this fork's extensions):

```{r}
if (!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("jaredsmurray/bcf")
```

On macOS without CRAN's gfortran toolchain, the link step fails with
`ld: library 'emutls_w' not found` because Makeconf points `FLIBS` at the
missing `/opt/gfortran` (the package itself has no Fortran). Install with
the override:

```sh
MAKEFLAGS='FLIBS=' R CMD INSTALL .
```
