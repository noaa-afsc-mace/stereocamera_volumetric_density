% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/detection.R
\name{fit_single_logistic_glm}
\alias{fit_single_logistic_glm}
\title{single logistic fitting - convenience function that rolls all three steps of fitting a function into one}
\usage{
fit_single_logistic_glm(
  target_ranges,
  vol_func,
  nbins,
  plotting = FALSE,
  exp_dens = NULL
)
}
\arguments{
\item{target_ranges}{vector of ranges for a given target}

\item{vol_func}{volume function as derived from the get_volume method}

\item{nbins}{number of bins to use for density histogram}

\item{plotting}{boolean flag to indicate whether to make plots}

\item{start_params}{Vector of parameters for single logistic (scale,r50,slope)}
}
\description{
single logistic fitting - convenience function that rolls all three steps of fitting a function into one
}
