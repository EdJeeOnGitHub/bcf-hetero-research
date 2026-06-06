\name{bcf_hetero}
\alias{bcf_hetero}
\title{Fit Bayesian Causal Forests With Shared Heteroscedastic Residual Variance}
\usage{
bcf_hetero(
  y,
  z,
  x_control,
  x_moderate = x_control,
  pihat,
  ...,
  x_variance = x_control,
  vartree = list(num_trees = 40, nu = 10, lambda = NULL, numcut = 100,
    sparse = FALSE),
  variance_model = "shared"
)
}
\arguments{
\item{y}{Response vector.}
\item{z}{Binary treatment vector.}
\item{x_control}{Design matrix for the prognostic function.}
\item{x_moderate}{Design matrix for the treatment effect function.}
\item{pihat}{Estimated propensity score.}
\item{...}{Additional arguments passed to \code{bcf()}.}
\item{x_variance}{Design matrix for the residual variance function.}
\item{vartree}{List of scalar product-of-trees variance prior parameters.}
\item{variance_model}{Currently only \code{"shared"} is implemented.}
}
\value{
A fitted \code{bcf_hetero} object. The object contains posterior draw matrices
\code{mu}, \code{tau}, \code{sigma2}, and \code{sigma}, where \code{sigma}
is \code{sqrt(sigma2)}.
}
\description{
Fits a BCF model with a shared covariate-dependent residual variance function
\eqn{\sigma_i^2 = d(X_i)}. The residual variance is sampled jointly with the
mean and treatment-effect trees, and BCF tree updates use observation-specific
residual precision.
}
\details{
When tree output is enabled with \code{no_output = FALSE}, serialized variance
trees are written alongside the existing BCF tree files. These can be used by
\code{predict(..., type = "sigma2")} to obtain out-of-sample residual variance
draws.
}
