\name{kde}

\alias{kde}

\title{
  Kernel Density Estimator over a Grid of Points
}

\description{
  Given a point cloud \code{X} (\eqn{n} points), this function computes the Kernel Density Estimator over a grid of points. The kernel is a Gaussian Kernel with smoothing parameter \code{h}. For each \eqn{x \in R^d}, the Kernel Density estimator is defined as
  \deqn{
    p_X(x)=\frac{1}{n (\sqrt{2 \pi} h )^d} \sum_{i=1}^n \exp \left( \frac{- \Vert x-X_i \Vert_2^2}{2h^2} \right).
  }   
}

\usage{
  kde(X, Grid, h, weight = 1, printProgress = FALSE)
}

\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates of points used in the kernel density estimation process, where \eqn{n} is the number of points and \eqn{d} is the dimension.
  }

  \item{Grid}{
    an \eqn{m} by \eqn{d} matrix of coordinates, where \eqn{m} is the number of points in the grid.
  }

  \item{h}{
    number: the smoothing paramter of the Gaussian Kernel.
  }

  \item{weight}{
    either a number, or a vector of length \eqn{n}. If it is a number, then same weight is applied to each points of \code{X}. If it is a vector, \code{weight} represents weights of each points of \code{X}.
  }

  \item{printProgress}{
    if TRUE a progress bar is printed. Default is \code{FALSE}.
  }
}

\value{
  \code{kde} returns a vector of length \eqn{m} (the number of points in the grid) containing the value of the kernel density estimator for each point in the grid.
}

\references{
  Larry Wasserman (2004), "All of statistics: a concise course in statistical inference", Springer.

  Brittany T. Fasy, Fabrizio Lecci, Alessandro Rinaldo, Larry Wasserman, Sivaraman Balakrishnan, and Aarti Singh. (2013), "Statistical Inference For Persistent Homology: Confidence Sets for Persistence Diagrams", (arXiv:1303.7117). To appear, Annals of Statistics.
}

\author{
  Jisu Kim, Fabrizio Lecci
}

\seealso{
  \code{\link{kernelDist}}, \code{\link{distFct}}, \code{\link{dtm}}
}

\examples{
## Generate Data from the unit circle
n <- 300
X <- circleUnif(n)

## Construct a grid of points over which we evaluate the function
by <- 0.065
Xseq <- seq(-1.6, 1.6, by=by)
Yseq <- seq(-1.7, 1.7, by=by)
Grid <- expand.grid(Xseq,Yseq)

## kernel density estimator
h <- 0.3
KDE <- kde(X, Grid, h)
}

\keyword{nonparametric}
