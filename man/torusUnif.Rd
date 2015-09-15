\name{torusUnif}

\alias{torusUnif}

\title{
  Uniform Sample From The 3D Torus
}

\description{
  This function samples \code{n} points from the 3D torus, uniformly with respect to its surface.
}

\usage{
  torusUnif(n, a, c)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{n}{
    an integer specifying the number of points in the sample.
  }

  \item{a}{
    the radius of the torus tube.
  }

  \item{c}{
    the radius from the center of the hole to the center of the torus tube.
  }
}

\details{
  This function is an implementation of Algorithm 1 in the reference.
}

\value{
  \code{torusUnif} returns an \code{n} by 3 matrix of coordinates.
}

\references{
  Persi Diaconis, Susan Holmes, and Mehrdad Shahshahani, (2013), "Sampling from a manifold." Advances in Modern Statistical Theory and Applications: A Festschrift in honor of Morris L. Eaton. Institute of Mathematical Statistics, 102-125.
}

\author{
  Fabrizio Lecci
}

\seealso{
  \code{\link{circleUnif}},\code{\link{sphereUnif}}
}

\examples{
X <- torusUnif(300, a = 1.8, c = 5)
plot(X)
}

\keyword{datagen}
