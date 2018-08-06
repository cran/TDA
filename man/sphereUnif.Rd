\name{sphereUnif}

\alias{sphereUnif}

\title{
  Uniform Sample From The Sphere \eqn{S^d}
}

\description{
  The function \code{sphereUnif} samples \code{n} points from the sphere \eqn{S^d} of radius \code{r} embedded in \eqn{R^{d+1}}, uniformly with respect to the volume measure of the sphere.
}

\usage{
sphereUnif(n, d, r = 1)
}

\arguments{
  \item{n}{
    an integer specifying the number of points in the sample.
  }

  \item{d}{
    an integer specifying the dimension of the sphere \eqn{S^d}
  }

  \item{r}{
    a numeric variable specifying the radius of the sphere. The default value is \code{1}.
  }
}

\value{
  The function \code{sphereUnif} returns an \code{n} by 2 matrix of coordinates.
}

\author{
  Jisu Kim
}

\note{
  When \code{d = 1}, this function is same as using \code{\link{circleUnif}}.
}

\seealso{
  \code{\link{circleUnif}}, \code{\link{torusUnif}}
}

\examples{
X <- sphereUnif(n = 100, d = 1, r = 1)
plot(X)
}

\keyword{datagen}
