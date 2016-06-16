\name{hausdInterval}

\alias{hausdInterval}

\title{
  Subsampling Confidence Interval for the Hausdorff Distance between a Manifold and a Sample
}

\description{
  \code{hausdInterval} computes a confidence interval for the Hausdorff distance between a point cloud \code{X} and the underlying manifold from which \code{X} was sampled. See Details and References.
}

\usage{
hausdInterval(X, m, B = 30, alpha = 0.05, parallel = FALSE,
              printProgress = FALSE)
}

\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates of sampled points.
  }

  \item{m}{
    the size of the subsamples.
  }

  \item{B}{
    the number of subsampling iterations. The default value is \code{30}.
  }

  \item{alpha}{
    \code{hausdInterval} returns a (\code{1-alpha}) confidence interval. The default value is \code{0.05}.
  }

  \item{parallel}{
    logical: if \code{TRUE}, the iterations are parallelized, using the library \code{parallel}. The default value is \code{FALSE}.
  }

  \item{printProgress}{
    if \code{TRUE}, a progress bar is printed. The default value is \code{FALSE}.
  }
}

\details{
  For \code{B} times, the subsampling algorithm subsamples \code{m} points of \code{X} (without replacement) and computes the Hausdorff distance between the original sample \code{X} and the subsample. The result is a sequence of \code{B} values. Let \eqn{q} be the (\code{1-alpha}) quantile of these values and let \eqn{c = 2 * q}. The interval \eqn{[0, c]} is a valid  (\code{1-alpha}) confidence interval for the Hausdorff distance between \code{X} and the underlying manifold, as proven in (Fasy, Lecci, Rinaldo, Wasserman, Balakrishnan, and Singh, 2013, Theorem 3).
}

\value{ 
  The function \code{hausdInterval} returns a number \eqn{c}. The confidence interval is \eqn{[0, c]}. 
}

\references{
  Fasy BT, Lecci F, Rinaldo A, Wasserman L, Balakrishnan S, Singh A (2013). "Statistical Inference For Persistent Homology: Confidence Sets for Persistence Diagrams." (arXiv:1303.7117). Annals of Statistics.
}

\author{
  Fabrizio Lecci
}

\seealso{
  \code{\link{bootstrapBand}}
}

\examples{
X <- circleUnif(1000)
interval <- hausdInterval(X, m = 800)
print(interval)
}

\keyword{nonparametric}
\keyword{htest}
