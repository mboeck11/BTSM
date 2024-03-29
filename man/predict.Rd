% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/predict.R
\name{predict}
\alias{predict}
\title{Predictions}
\usage{
predict(
  object,
  ...,
  n.ahead = 4,
  quantiles = c(0.05, 0.1, 0.16, 0.5, 0.84, 0.9, 0.95),
  applyfun = NULL,
  cores = NULL,
  save.store = FALSE,
  verbose = TRUE
)
}
\arguments{
\item{object}{an object of class \code{bvar}.}

\item{...}{additional arguments.}

\item{n.ahead}{the forecast horizon.}

\item{quantiles}{posterior quantiles to be computed.}

\item{applyfun}{parallelization}

\item{cores}{number of cores}

\item{save.store}{If set to \code{TRUE} the full distribution is returned. Default is set to \code{FALSE} in order to save storage.}

\item{verbose}{If set to \code{FALSE} it suppresses printing messages to the console.}
}
\description{
Predictions
}
