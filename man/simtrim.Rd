\name{simtrim}
\alias{simtrim}

\title{
Simulation and Statistics for Trimmed Mean Estimation of Weibull Distribution
}

\description{
This function conducts simulations and calculates various statistics related to the trimmed mean estimation for a Weibull distribution.
}

\usage{
simtrim(n, shape, scale, p, S)
}

\arguments{
  \item{n}{
    The sample size (integer greater than 1) for each simulation.
  }

  \item{shape}{
    The shape parameter (positive numeric) of the Weibull distribution.
  }

  \item{scale}{
    The scale parameter (positive numeric) of the Weibull distribution.
  }

  \item{p}{
    The proportion (numeric) of values to trim from the upper end of the distribution for the trimmed mean estimation. Value should be between 0 and 1.
  }

  \item{S}{
    The number of simulations (positive integer) to conduct.
  }
}

\details{
The function calculates the the trimmed mean (\code{mu_p}) for \code{S} simulations. On each simulation, \code{n} observations are generated from a Weibull distribution with parameters \code{shape} and \code{scale}. The output is the mean, bias, variance, and MSE of the standard untrimmed sample mean (\code{mu_0}, computed using properties of the Weibull distribution) and the trimmed mean (\code{mu_p}, estimated based on the simulations).
}

\value{
A list containing the following components:

\item{mu_0}{
  The mean of the untrimmed sample average (based on Weibull distribution)
}

\item{bias_mu_0}{
  The bias of the untrimmed sample average (always 0)
}

\item{var_mu_0}{
  The variance of the untrimmed sample average (based on Weibull distribution)
}

\item{MSE_mu_0}{
  The MSE of the untrimmed sample average (based on Weibull distribution)
}


\item{mu_p}{
  The mean of the trimmed sample average (based on simulations)
}

\item{bias_mu_p}{
  The bias of the trimmed sample average (based on simulations)
}

\item{var_mu_p}{
  The variance of the trimmed sample average (based on simulations)
}

\item{MSE_mu_p}{
  The MSE of the trimmed sample average (based on simulations)
}



}
