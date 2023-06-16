\name{wei_mu_mle}
\alias{wei_mu_mle}

\title{
Maximum Likelihood Estimate of Mean for Weibull Distribution
}

\description{
This function calculates the maximum likelihood estimate (MLE) for the mean of a Weibull distribution using the negative log-likelihood function.
}

\usage{
wei_mu_mle(y)
}

\arguments{
  \item{y}{
    A numeric vector of observed data.
  }
}

\details{
The function estimates the mean of a Weibull distribution using the maximum likelihood estimation method. It utilizes the negative log-likelihood function and the \code{nlminb} optimization function from the \code{stats} package to find the parameter values that maximize the likelihood of the observed data.
}

\value{
The maximum likelihood estimate (MLE) of the mean for the Weibull distribution.
}

\examples{
# Example usage
y <- c(1.2, 1.8, 2.5) # Example observed data
wei_mu_mle(y)
}
