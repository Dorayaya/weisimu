\name{simtrim_by}
\alias{simtrim_by}

\title{
Plot Simulations for Trimmed Mean Estimation of Weibull Distribution
}

\description{
This function conducts simulations and plots MSE for trimmed vs untrimmed means.
}

\usage{
simtrim_by(n, shape, scale, p, S, plot=TRUE, lty=c(2,1), lwd=c(2,2), col=c("black","black"), leg.pos="topright", ...)
}

\arguments{
  \item{n}{
    See details.
  }

  \item{shape}{
    See details.
  }

  \item{scale}{
    See details.
  }

  \item{p}{
    See details.
  }

  \item{S}{
    The number of simulations (positive integer) to conduct.
  }

  \item{plot}{
    If TRUE (default) plot the results.
  }

  \item{lty}{
    If plot generated, this controls line type for trimmed and untrimmed groups. Specified as numeric vector of length 2. See base R plotting for details.
  }

  \item{lwd}{
    If plot generated, this controls line width for trimmed and untrimmed groups. Specified as numeric vector of length 2. See base R plotting for details.
  }

  \item{col}{
    If plot generated, this controls line colour for trimmed and untrimmed groups. Specified as character vector of length 2. See base R plotting for details.
  }

  \item{leg.pos}{
    If plot generated, this is a string denoting where the legend should be placed. See \code{?legend} for details.
  }

  \item{...}{
    Any additional parameters are sent to the \code{plot} function.
  }

}

\details{
\code{n}, \code{shape}, \code{scale}, \code{p} are as described in \code{simtrim}. However, exactly one must be specified as a vector. For each value in the vector, we run a simulation as in \code{simtrim}. The parameter which is a vector will be used as the x-axis of the plot, while the y-axis is the MSE.
}

\value{
The same output as \code{simtrim}, but repeated for each value of the parameter which is specified as a vector.
}
