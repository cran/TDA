\name{knnDE}

\alias{knnDE}

\title{
  k Nearest Neighbors Density Estimator over a Grid of Points
}

\description{
  Given a point cloud \code{X} (\eqn{n} points), this function computes the k Nearest Neighbors Density Estimator over a grid of points. For each \eqn{x \in R^d}, the knn Density Estimator is defined by
  \deqn{
    p_X(x)=\frac{k}{n \; v_d \; r_k^d(x)},
  }   
  where \eqn{v_n} is the volume of the Euclidean \eqn{d} dimensional unit ball and \eqn{r_k^d(x)} is the Euclidean distance form point x to its \eqn{k}'th closest neighbor.
}

\usage{
knnDE(X, Grid, k)
}

\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates of points used in the density estimation process, where \eqn{n} is the number of points and \eqn{d} is the dimension.
  }

  \item{Grid}{
    an \eqn{m} by \eqn{d} matrix of coordinates, where \eqn{m} is the number of points in the grid.
  }

  \item{k}{
    number: the smoothing paramter of the k Nearest Neighbors Density Estimator.
  }
}

\value{
  \code{knnDE} returns a vector of length \eqn{m} (the number of points in the grid) containing the value of the knn Density Estimator for each point in the grid.
}

% \references{
% Larry Wasserman (2004), "All of statistics: a concise course in statistical inference", Springer.
% }

\author{
  Fabrizio Lecci
}

\seealso{
  \code{\link{kde}},\code{\link{kernelDist}}, \code{\link{distFct}}, \code{\link{dtm}}
}

\examples{
## Generate Data from the unit circle
n <- 300
X <- circleUnif(n)

## Construct a grid of points over which we evaluate the function
by <- 0.065
Xseq <- seq(-1.6, 1.6, by = by)
Yseq <- seq(-1.7, 1.7, by = by)
Grid <- expand.grid(Xseq, Yseq)

## kernel density estimator
k <- 50
KNN <- knnDE(X, Grid, k)
}

\keyword{nonparametric}

