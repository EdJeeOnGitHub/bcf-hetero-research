\name{bcf_hetero_diagnostics}
\alias{bcf_hetero_diagnostics}
\title{Diagnostics for Heteroscedastic BCF Simulation Comparisons}
\usage{
bcf_hetero_diagnostics(
  fit,
  fit_reference = NULL,
  true_sigma2 = NULL,
  alpha = NULL,
  top_k = NULL,
  exponent_bounds = c(-700, 700)
)
}
\arguments{
\item{fit}{A fitted object returned by \code{bcf_hetero()}.}
\item{fit_reference}{Optional fitted standard BCF object used for overlap and rank comparisons.}
\item{true_sigma2}{Optional vector of true residual variances.}
\item{alpha}{Optional CARA risk-aversion parameter. If supplied, variance-aware scores are computed.}
\item{top_k}{Optional number of observations selected for overlap/Jaccard comparisons.}
\item{exponent_bounds}{Numeric length-two bounds used by \code{bcf_hetero_cara_score()}.}
}
\value{
A list containing posterior summaries and any requested diagnostics, including
\code{sigma2_correlation}, \code{cara_score}, \code{rank_correlation},
\code{rank_change}, and top-k overlap/Jaccard fields when the corresponding
inputs are supplied.
}
\description{
Computes the heteroscedastic BCF diagnostics used by simulation comparisons:
posterior mean residual variance, optional correlation with true residual
variance, variance-aware CARA scores, rank changes relative to a reference BCF
fit, and top-k overlap/Jaccard.
}
