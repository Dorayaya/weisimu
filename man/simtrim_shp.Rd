\name{simtrim_shp}
\alias{simtrim_shp}
\title{Simulate and Calculate Mean Squared Errors (MSE) for Various Shape Parameters}
\description{
  Simulate data from a Weibull distribution with different shape parameters and calculate the MSEs for mu_0 and mu_p estimators.
}
\usage{
  simtrim_shp(n, shape = NA, scale, p, S)
}
\arguments{
  \item{n}{An integer specifying the sample size.}
  \item{shape}{An optional numeric vector representing the range of shape parameters. If not provided, the function will use the default range from 0.05 to 1 with an increment of 0.01.}
  \item{scale}{A numeric value representing the scale parameter of the Weibull distribution.}
  \item{p}{A numeric value indicating the trimming proportion.}
  \item{S}{An integer specifying the number of simulations to be conducted.}
}
\value{
  A list containing the shape parameters, MSE values for mu_0, and MSE values for mu_p.
}
\examples{
  # Simulate and calculate MSEs for different shape parameters
  results <- simtrim_shp(n = 100, scale = 3, p = 0.1, S = 100)
  
  # Print the shape parameters
  print(results[[1]])
  
  # Print the MSEs for mu_0
  print(results[[2]])
  
  # Print the MSEs for mu_p
  print(results[[3]])
}
