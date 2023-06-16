\name{mean}
\alias{mean}

\title{
Trimmed Mean Function with Upper Trimming
}

\description{
This function calculates the mean of a numeric vector, with the option to trim values above a certain threshold.
}

\usage{
mean(x, trim = 0, na.rm = FALSE, trim.upper = FALSE, ...)
}

\arguments{
  \item{x}{
    A numeric vector.
  }

  \item{trim}{
    A numeric value specifying the proportion of values to trim from both ends of the vector. Default is 0, meaning no trimming is performed.
  }

  \item{na.rm}{
    A logical value indicating whether missing values should be removed from the vector before calculating the mean. Default is \code{FALSE}.
  }

  \item{trim.upper}{
    A logical value indicating whether to trim values above the upper quantile threshold. If \code{TRUE}, values above the quantile threshold will be excluded from the mean calculation. If \code{FALSE}, the regular mean calculation will be performed. Default is \code{FALSE}.
  }

  \item{...}{
    Additional arguments to be passed to the base \code{mean} function.
  }
}

\details{
The function calculates the mean of the input vector. If \code{trim.upper} is \code{FALSE}, the mean is calculated without any trimming using the base \code{mean} function. If \code{trim.upper} is \code{TRUE}, values above the upper quantile threshold are trimmed before calculating the mean.
}

\value{
The calculated mean value.
}

\examples{
# Example usage
x <- c(1, 2, 3, 4, 5)
mean(x) # Regular mean
mean(x, trim = 0.2) # Trim 20% of values from both ends
mean(x, trim.upper = TRUE) # Trim values above the upper quantile threshold
mean(x, trim = 0.2, trim.upper = TRUE) # Trim values from both ends and above the upper quantile threshold
}
