\name{simtrim_ss}
\alias{simtrim_ss}
\title{Simulate and Calculate Mean Squared Errors (MSE) for Various Sample Sizes}
\description{
  Simulate data from a Weibull distribution with different sample sizes and calculate the MSEs for mu_0 and mu_p estimators.
}
\usage{
  simtrim_ss(n = NA, shape, scale, p, S)
}
\arguments{
  \item{n}{An optional integer specifying the sample size. If not provided, the function will use the default range from 20 to 2000.}
  \item{shape}{A numeric value representing the shape parameter of the Weibull distribution.}
  \item{scale}{A numeric value representing the scale parameter of the Weibull distribution.}
  \item{p}{A numeric value indicating the trimming proportion.}
  \item{S}{An integer specifying the number of simulations to be conducted.}
}
\value{
  A list containing the sample sizes, MSE values for mu_0, and MSE values for mu_p.
}
\examples{
  # Simulate and calculate MSEs for different sample sizes
  results <- simtrim_ss(shape = 2, scale = 3, p = 0.1, S = 100)
  
  # Print the sample sizes
  print(results[[1]])
  
  # Print the MSEs for mu_0
  print(results[[2]])
  
  # Print the MSEs for mu_p
  print(results[[3]])
}
