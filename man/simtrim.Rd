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
    The sample size (integer) for each simulation.
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
    The number of simulations (integer) to conduct.
  }
}

\details{
The function conducts \code{S} simulations of a Weibull distribution with parameters \code{shape} and \code{scale}. For each simulation, it calculates the regular mean (\code{mu_0}), the trimmed mean (\code{mu_p}), and the mean estimated using maximum likelihood estimation (\code{mu_mle}). It also calculates various statistics such as bias, variance, and mean squared error (MSE) for each estimation method.
}

\value{
A list containing the following components:

\item{true_mean}{
  The true mean of the Weibull distribution.
}

\item{bias}{
  A list of bias values for each estimation method:
  \begin{itemize}
    \item \code{mu_0}: Bias for the regular mean estimation.
    \item \code{mu_p}: Bias for the trimmed mean estimation.
    \item \code{mu_mle}: Bias for the mean estimation using maximum likelihood estimation.
  \end{itemize}
}

\item{vars}{
  A list of variance values for each estimation method:
  \begin{itemize}
    \item \code{mu_0}: Variance for the regular mean estimation.
    \item \code{mu_p}: Variance for the trimmed mean estimation.
    \item \code{mu_mle}: Variance for the mean estimation using maximum likelihood estimation.
  \end{itemize}
}

\item{MSE}{
  A list of mean squared error (MSE) values for each estimation method:
  \begin{itemize}
    \item \code{mu_0}: MSE for the regular mean estimation.
    \item \code{mu_p}: MSE for the trimmed mean estimation.
    \item \code{mu_mle}: MSE for the mean estimation using maximum likelihood estimation.
  \end
