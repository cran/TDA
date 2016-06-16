\name{alphaShapeDiag}

\alias{alphaShapeDiag}

\title{
  Persistence Diagram of Alpha Shape in 3d
}

\description{
  The function \code{alphaShapeDiag} computes the persistence diagram of the alpha shape filtration built on top of a point cloud in 3 dimension.
}

\usage{
alphaShapeDiag(X, library = "GUDHI", printProgress = FALSE)
}

\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates, used by the function \code{FUN}, where \eqn{n} is the number of points stored in \code{X} and \eqn{d} is the dimension of the space.
    Currently \eqn{d} should be 3.
  }

  \item{library}{
    Currently only the library "GUDHI" can be used to compute the persistence diagram.
  }

  \item{printProgress}{
    if \code{TRUE}, a progress bar is printed. The default value is \code{FALSE}.
  }
}

\details{
  The function \code{alphaShapeDiag} constructs the alpha shape filtration, using the C++ library \href{ https://project.inria.fr/gudhi/software }{GUDHI}.
  See refereneces.
}

\value{
  The function \code{alphaShapeDiag} returns a list with the following elements:
  \item{diagram}{
    an object of class \code{diagram}, a \eqn{P} by 3 matrix, where \eqn{P} is the number of points in the resulting persistence diagram. The first column stores the dimension of each feature (0 for components, 1 for loops, 2 for voids, etc). Second and third columns are Birth and Death of the features.
  }
}

\references{
  Fischer K (2005). "Introduction to Alpha Shapes."

  Edelsbrunner H, Mucke EP (1994). "Three-dimensional Alpha Shapes." ACM Trans. Graph.

  Maria C (2014). "GUDHI, Simplicial Complexes and Persistent Homology Packages." \url{ https://project.inria.fr/gudhi/software/ }

  Morozov D (2008). "Homological Illusions of Persistence and Stability." \url{ http://www.mrzv.org/publications/thesis/phd/ }

  Edelsbrunner H, Harer J (2010). "Computational topology: an introduction." American Mathematical Society.
}

\author{
  Jisu Kim and Vincent Rouvreau
}

\seealso{
  \code{\link{summary.diagram}}, \code{\link{plot.diagram}}, \code{\link{alphaComplexDiag}}, \code{\link{gridDiag}}, \code{\link{ripsDiag}}
}

\examples{
# input data generated from sphere
XX <- sphereUnif(n = 500, d = 2)

# persistence diagram of alpha shape
DiagAlphaShape <- alphaShapeDiag(X = XX, printProgress = TRUE)

# plot
plot(DiagAlphaShape[["diagram"]])
}

\keyword{methods}
