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

Install the latest development version from GitHub:

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
