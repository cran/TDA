\name{bootstrapBand}

\alias{bootstrapBand}

\title{
  Bootstrap Confidence Band
}

\description{
  The function \code{bootstrapBand} computes a uniform symmetric confidence band around a function of the data \code{X}, evaluated on a \code{Grid}, using the bootstrap algorithm. See Details and References.
}

\usage{
bootstrapBand(
    X, FUN, Grid, B = 30, alpha = 0.05, parallel = FALSE,
    printProgress = FALSE, weight = NULL, ...)
}

\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates of points used by the function \code{FUN}, where \eqn{n} is the number of points and \eqn{d} is the dimension.
  }

  \item{FUN}{
    a function whose inputs are an \eqn{n} by \eqn{d} matrix of coordinates \code{X}, an \eqn{m} by \eqn{d} matrix of coordinates \code{Grid} and returns a numeric vector of length \eqn{m}. For example see \code{\link{distFct}}, \code{\link{kde}}, and \code{\link{dtm}}  which compute the distance function, the kernel density estimator and the distance to measure over a grid of points, using the input \code{X}.
  }

  \item{Grid}{
    an \eqn{m} by \eqn{d} matrix of coordinates, where \eqn{m} is the number of points in the grid, at which \code{FUN} is evaluated.
  }

  \item{B}{
    the number of bootstrap iterations.
  }

  \item{alpha}{
    \code{bootstrapBand} returns a (\code{1-alpha}) confidence band. The default value is \code{0.05}.
  }

  \item{parallel}{
    logical: if \code{TRUE} the bootstrap iterations are parallelized, using the library \code{parallel}. The default value is \code{FALSE}.
  }

  \item{printProgress}{
    if \code{TRUE}, a progress bar is printed. The default value is \code{FALSE}.
  }

  \item{weight}{
    either NULL, a number, or a vector of length \eqn{n}. If it is NULL, weight is not used. If it is a number, then same weight is applied to each points of \code{X}. If it is a vector, \code{weight} represents weights of each points of \code{X}. The default value is \code{NULL}.
  }

  \item{...}{
    additional parameters for the function \code{FUN}.
  }
}

\details{
  First, the input function \code{FUN} is evaluated on the \code{Grid} using the original data \code{X}. Then, for \code{B} times, the bootstrap algorithm subsamples \code{n} points of \code{X} (with replacement), evaluates the function \code{FUN} on the \code{Grid} using the subsample, and computes the \eqn{\ell_\infty}{l_\infty} distance between the original function and the bootstrapped one. The result is a sequence of \code{B} values. The (\code{1-alpha}) confidence band is constructed by taking the (\code{1-alpha}) quantile of these values.
}

\value{
  The function \code{bootstrapBand} returns a list with the following elements:

  \item{width}{
    number: (\code{1-alpha}) quantile of the values computed by the bootstrap algorithm. It corresponds to half of the width of the unfiorm confidence band; that is, \code{width} is the distance of the upper and lower limits of the band from the function evaluated using the original dataset \code{X}.
  }

  \item{fun}{
    a numeric vector of length \eqn{m}, storing the values of the input function \code{FUN}, evaluated on the \code{Grid} using the original data \code{X}.
  }

  \item{band}{
    an \eqn{m} by 2 matrix that stores the values of the lower limit of the confidence band (first column) and upper limit of the confidence band (second column), evaluated over the \code{Grid}.
  }
}

\references{
  Wasserman L (2004). "All of statistics: a concise course in statistical inference." Springer.

  Fasy BT, Lecci F, Rinaldo A, Wasserman L, Balakrishnan S, Singh A (2013). "Statistical Inference For Persistent Homology: Confidence Sets for Persistence Diagrams." (arXiv:1303.7117). Annals of Statistics.

  Chazal F, Fasy BT, Lecci F, Michel B, Rinaldo A, Wasserman L (2014). "Robust Topological Inference: Distance-To-a-Measure and Kernel Distance." Technical Report. 
}

\author{
  Jisu Kim and Fabrizio Lecci
}

\seealso{
  \code{\link{kde}}, \code{\link{dtm}}
}

\examples{
# Generate data from mixture of 2 normals.
n <- 2000
X <- c(rnorm(n / 2), rnorm(n / 2, mean = 3, sd = 1.2))

# Construct a grid of points over which we evaluate the function
by <- 0.02
Grid <- seq(-3, 6, by = by)

## bandwidth for kernel density estimator
h <- 0.3
## Bootstrap confidence band
band <- bootstrapBand(X, kde, Grid, B = 80, parallel = FALSE, alpha = 0.05,
                      h = h)

plot(Grid, band[["fun"]], type = "l", lwd = 2,
     ylim = c(0, max(band[["band"]])), main = "kde with 0.95 confidence band")
lines(Grid, pmax(band[["band"]][, 1], 0), col = 2, lwd = 2)
lines(Grid, band[["band"]][, 2], col = 2, lwd = 2)
}

\keyword{nonparametric}
\keyword{htest}
