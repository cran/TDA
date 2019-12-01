\name{alphaShapeFiltration}

\alias{alphaShapeFiltration}

\title{
  Alpha Shape Filtration in 3d
}

\description{
  The function \code{alphaShapeFiltration} computes the alpha shape filtration built on top of a point cloud in 3 dimension.
}

\usage{
alphaShapeFiltration(
    X, library = "GUDHI", printProgress = FALSE)
}

\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates, used by the function \code{FUN}, where \eqn{n} is the number of points stored in \code{X} and \eqn{d} is the dimension of the space.
    Currently \eqn{d} should be 3.
  }

  \item{library}{
    a string specifying which library to compute the Alpha Shape filtration. The user can use the library \code{"GUDHI"}, and is also the default value.
  }

  \item{printProgress}{
    if \code{TRUE}, a progress bar is printed. The default value is \code{FALSE}.
  }
}

\details{
  The function \code{alphaShapeFiltration} constructs the alpha shape filtration, using the C++ library \href{ https://project.inria.fr/gudhi/software }{GUDHI}.
  See refereneces.
}

\value{
  The function \code{alphaShapeFiltration} returns a list with the following elements:

  \item{cmplx}{
    a list representing the complex. Its i-th element represents the vertices of i-th simplex.
  }

  \item{values}{
    a vector representing the filtration values. Its i-th element represents the filtration value of i-th simplex.
  }

  \item{increasing}{
    a logical variable indicating if the filtration values are in increasing order (\code{TRUE}) or in decreasing order (\code{FALSE}).
  }

  \item{coordinates}{
    a matrix representing the coordinates of vertices. Its i-th row represents the coordinate of i-th vertex.
  }
}

\references{
  Fischer K (2005). "Introduction to Alpha Shapes."

  Edelsbrunner H, Mucke EP (1994). "Three-dimensional Alpha Shapes." ACM Trans. Graph.

  Maria C (2014). "GUDHI, Simplicial Complexes and Persistent Homology Packages." \url{ https://project.inria.fr/gudhi/software/ }

  Morozov D (2008). "Homological Illusions of Persistence and Stability."

  Edelsbrunner H, Harer J (2010). "Computational topology: an introduction." American Mathematical Society.
}

\author{
  Jisu Kim and Vincent Rouvreau
}

\seealso{
  \code{\link{alphaShapeDiag}}, \code{\link{filtrationDiag}}
}

\examples{
# input data generated from sphere
X <- sphereUnif(n = 20, d = 2)

# alpha shape filtration
FltAlphaShape <- alphaShapeFiltration(X = X, printProgress = TRUE)
}

\keyword{methods}
