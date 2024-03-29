% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NBFitPlot.R
\name{NBFitPlot}
\alias{NBFitPlot}
\title{Plot the goodness of fitting}
\usage{
NBFitPlot(dat, seed = 1)
}
\arguments{
\item{dat}{Expression counts vector for a single cell.}

\item{seed}{Random seed for the simulated curve from theoretical models.}
}
\value{
A plot with observed and theoretical distribution of the expression counts.
}
\description{
This function fit the Negative Binomial model on one single cell and plot the
frequecy of oberseved counts and the theoretical probability from the NB model.
This function helps evaluate the goodness of fitting of the model to the dataset
before applying.
}
