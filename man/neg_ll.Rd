\name{neg_ll}
\alias{neg_ll}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Log-Likelihood Function for Weibull Distribution
}
\description{
This function calculates the negative log-likelihood for the Weibull distribution.
}
\usage{
neg_ll(theta, y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
A numeric vector containing the shape and scale parameters of the Weibull distribution.
}
}
\iten{y}{
 A numeric vector of observed data.
}
\details{The negative log-likelihood is computed using the formula for the Weibull distribution. It evaluates the goodness of fit of the Weibull distribution to the observed data.
}}
\value{
A numeric value representing the negative log-likelihood.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Gradon Nicholls, Xiaoya Wang
}
\note{
This function assumes that the observed data `y` follows a Weibull distribution.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{dweibull}}, \code{\link{pweibull}}, \code{\link{qweibull}}, \code{\link{rweibull}}
}

\examples{
# Example usage
theta <- c(2, 1) # Example values for shape and scale parameters
y <- c(1.2, 1.8, 2.5) # Example observed data
neg_ll(theta, y)
}
