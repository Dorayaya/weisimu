\name{simtrim_p}
\alias{simtrim_p}
\title{Simulate and Calculate Mean Squared Errors (MSE) for Various Trimming Proportions}
\description{
  Simulate data from a Weibull distribution with different trimming proportions and calculate the MSEs for mu_0 and mu_p estimators.
}
\usage{
  simtrim_p(n, shape, scale, p = NA, S)
}
\arguments{
  \item{n}{An integer specifying the sample size.}
  \item{shape}{A numeric value representing the shape parameter of the Weibull distribution.}
  \item{scale}{A numeric value representing the scale parameter of the Weibull distribution.}
  \item{p}{An optional numeric vector indicating the range of trimming proportions. If not provided, the function will use the default range from 0.01 to 0.20 with an increment of 0.01.}
  \item{S}{An integer specifying the number of simulations to be conducted.}
}
\value{
  A list containing the trimming proportions, MSE values for mu_0, and MSE values for mu_p.
}
\examples{
  # Simulate and calculate MSEs for different trimming proportions
  results <- simtrim_p(n = 100, shape = 2, scale = 3, S = 100)
  
  # Print the trimming proportions
  print(results[[1]])
  
  # Print the MSEs for mu_0
  print(results[[2]])
  
  # Print the MSEs for mu_p
  print(results[[3]])
}
