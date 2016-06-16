\name{kernelDist}

\alias{kernelDist}

\title{
  Kernel distance over a Grid of Points
}

\description{
  Given a point cloud \code{X}, the function \code{kernelDist} computes the kernel distance over a grid of points. The kernel is a Gaussian Kernel with smoothing parameter h:
  \deqn{
    K_h(x,y)=\exp\left( \frac{- \Vert x-y \Vert_2^2}{2h^2} \right).
  }{
    K_h (x,y) = exp( -(||x-y||^2)/(2h^2) ).
  }
  For each \eqn{x \in R^d}{x in R^d}, the Kernel distance is defined by
  \deqn{
    \kappa_X(x)=\sqrt{ \frac{1}{n^2} \sum_{i=1}^n\sum_{j=1}^n K_h(X_i, X_j) + K_h(x,x) - 2 \frac{1}{n} \sum_{i=1}^n K_h(x,X_i)  }. 
  }{
    \kappa_X (x) = \sqrt( 1/(n^2) \sum_{i=1}^n\sum_{j=1}^n K_h(X_i, X_j) + K_h(x,x) - 2/n \sum_{i=1}^n K_h(x,X_i) ). 
  }
}

\usage{
kernelDist(X, Grid, h, weight = 1, printProgress = FALSE)
}

\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates of points, where \eqn{n} is the number of points and \eqn{d} is the dimension.
  }

  \item{Grid}{
    an \eqn{m} by \eqn{d} matrix of coordinates, where \eqn{m} is the number of points in the grid.
  }

  \item{h}{
    number: the smoothing paramter of the Gaussian Kernel.
  }

  \item{weight}{
    either a number, or a vector of length \eqn{n}. If it is a number, then same weight is applied to each points of \code{X}. If it is a vector, \code{weight} represents weights of each points of \code{X}. The default value is \code{1}.
  }

  \item{printProgress}{
    if \code{TRUE}, a progress bar is printed. The default value is \code{FALSE}.
  }
}

\value{
  The function \code{kernelDist} returns a vector of lenght \eqn{m} (the number of points in the grid) containing the value of the Kernel distance for each point in the grid.
}

\references{
  Phillips JM, Wang B, Zheng Y (2013). "Geometric Inference on Kernel Density Estimates." arXiv:1307.7760.

  Chazal F, Fasy BT, Lecci F, Michel B, Rinaldo A, Wasserman L (2014). "Robust Topological Inference: Distance-To-a-Measure and Kernel Distance." Technical Report. 
}

\author{
  Jisu Kim and Fabrizio Lecci
}

\seealso{
  \code{\link{kde}}, \code{\link{dtm}}, \code{\link{distFct}}
}

\examples{
## Generate Data from the unit circle
n <- 300
X <- circleUnif(n)

## Construct a grid of points over which we evaluate the functions
by <- 0.065
Xseq <- seq(-1.6, 1.6, by = by)
Yseq <- seq(-1.7, 1.7, by = by)
Grid <- expand.grid(Xseq, Yseq)

## kernel distance estimator
h <- 0.3
Kdist <- kernelDist(X, Grid, h)
}

\keyword{nonparametric}