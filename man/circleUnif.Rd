\name{circleUnif}

\alias{circleUnif}

\title{
  Uniform Sample From The Circle
}

\description{
  The function \code{circleUnif} samples \code{n} points from the circle of radius \code{r}, uniformly with respect to the circumference length.
}

\usage{
  circleUnif(n, r = 1)
}

\arguments{
  \item{n}{
    an integer specifying the number of points in the sample.
  }

  \item{r}{
    a numeric variable specifying the radius of the circle. The default value is \code{1}.
  }
}

\value{
  \code{circleUnif} returns an \code{n} by 2 matrix of coordinates.
}

\author{
  Fabrizio Lecci
}

\note{
  Uniform sample from sphere of arbitrary dimension can be generated using \code{\link{sphereUnif}}.
}

\seealso{
  \code{\link{sphereUnif}}, \code{\link{torusUnif}}
}

\examples{
X <- circleUnif(100)
plot(X)
}

\keyword{datagen}
