\name{plotmse}
\alias{plotmse}
\title{Plot Mean Squared Errors (MSE)}
\description{
  Plot the Mean Squared Errors (MSE) for two estimators, mu_0 and mu_p, against a variable specified by type (sample size, shape, or trimming proportion).
}
\usage{
  plotmse(x, y1, y2, type)
}
\arguments{
  \item{x}{A numeric vector representing the variable values to be plotted.}
  \item{y1}{A numeric vector containing the MSE values for the mu_0 estimator.}
  \item{y2}{A numeric vector containing the MSE values for the mu_p estimator.}
  \item{type}{An integer indicating the type of variable to be plotted: 1 for sample size, 2 for shape, or 3 for trimming proportion.}
}
\examples{
  # Generate example data
  x <- 1:10
  y1 <- runif(10)
  y2 <- runif(10)
  
  # Plot MSE
  plotmse(x, y1, y2, type = 1)
}
