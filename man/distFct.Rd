\name{distFct}

\alias{distFct}

\title{
  Distance function
}

\description{
  The function \code{distFct} computes the distance between each point of a set \code{Grid} and the corresponding closest point of another set \code{X}.
}

\usage{
distFct(X, Grid)
}

\arguments{
  \item{X}{
    a numeric \eqn{m} by \eqn{d} matrix of coordinates in the space, where \eqn{m} is the number of points in \code{X} and \eqn{d} is the dimension of the space. 
    X is the set of points whose distance is being measured from a base grid.
  }

  \item{Grid}{
    a numeric \eqn{n} by \eqn{d} matrix of coordinates in the space, where \eqn{n} is the number of points in \code{Grid} and \eqn{d} is the dimension of the space. 
    Grid is the base set from which each point is compared to the closest point in X.
  }
}

\details{
  Given a set of points \code{X}, the distance function computed at \eqn{g} is defined as
  \deqn{
    d(g) = \inf_{x \in X} \| x-g \|_2
  }{
    d(g) = inf_{x in X} ||x-g||_2
  }
}

\value{
  The function \code{distFct} returns a numeric vector of length \eqn{n}, where \eqn{n} is the number of points stored in \code{Grid}. 
  Each value in V corresponds to the distance between a point in G and the nearest point in X.
}

\author{
  Fabrizio Lecci
}

\seealso{
  \code{\link{kde}},\code{\link{kernelDist}}, \code{\link{dtm}}
}

\examples{
## Generate Data from the unit circle
n <- 300
X <- circleUnif(n)

## Construct a grid of points over which we evaluate the function
interval <- 0.065
Xseq <- seq(-1.6, 1.6, by = interval)
Yseq <- seq(-1.7, 1.7, by = interval)
Grid <- expand.grid(Xseq, Yseq)

## distance fct
distance <- distFct(X, Grid)
}

\keyword{nonparametric}
