\name{bcf_hetero_cara_score}
\alias{bcf_hetero_cara_score}
\title{Variance-Aware CARA Treatment Scores From Heteroscedastic BCF Draws}
\usage{
bcf_hetero_cara_score(fit, alpha, exponent_bounds = c(-700, 700))
}
\arguments{
\item{fit}{A fitted object returned by \code{bcf_hetero()}.}
\item{alpha}{Nonnegative CARA risk-aversion parameter.}
\item{exponent_bounds}{Numeric length-two bounds applied before exponentiation
to avoid overflow.}
}
\value{
A vector of posterior mean treatment scores, one per observation.
}
\description{
Computes
\deqn{
E_s[
U(\mu_i^s + \tau_i^s, \sigma_i^{2,s}; \alpha)
-
U(\mu_i^s, \sigma_i^{2,s}; \alpha)
]
}
where
\deqn{
U(m, \sigma^2; \alpha) =
-\exp(-\alpha m + 0.5 \alpha^2 \sigma^2).
}
The exponent is clamped to \code{exponent_bounds} before exponentiation.
}
