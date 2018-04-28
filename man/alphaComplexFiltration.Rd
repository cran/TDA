\name{alphaComplexFiltration}

\alias{alphaComplexFiltration}

\title{
  Alpha Complex Filtration
}

\description{
  The function \code{alphaComplexFiltration} computes the alpha complex filtration built on top of a point cloud.
}

\usage{
alphaComplexFiltration(
    X, library = "GUDHI", printProgress = FALSE)
}

\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates, used by the function \code{FUN}, where \eqn{n} is the number of points stored in \code{X} and \eqn{d} is the dimension of the space.
  }

  \item{library}{
    a string specifying which library to compute the Alpha Complex filtration. The user can use the library \code{"GUDHI"}, and is also the default value.
  }

  \item{printProgress}{
    if \code{TRUE}, a progress bar is printed. The default value is \code{FALSE}.
  }
}

\details{
  The function \code{alphaComplexFiltration} constructs the alpha complex filtration, using the C++ library \href{ https://project.inria.fr/gudhi/software }{GUDHI}.
  See refereneces.
}

\value{
  The function \code{alphaComplexFiltration} returns a list with the following elements:

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

  Edelsbrunner H, Harer J (2010). "Computational topology: an introduction." American Mathematical Society.

  Rouvreau V (2015). "Alpha complex." In GUDHI User and Reference Manual. GUDHI Editorial Board. \url{ http://gudhi.gforge.inria.fr/doc/latest/group__alpha__complex.html }

  Edelsbrunner H, Kirkpatrick G, Seidel R (1983). "On the shape of a set
of points in the plane." IEEE Trans. Inform. Theory.

  Maria C (2014). "GUDHI, Simplicial Complexes and Persistent Homology Packages." \url{ https://project.inria.fr/gudhi/software/ }

}

\author{
  Jisu Kim and Vincent Rouvreau
}

\seealso{
  \code{\link{alphaComplexDiag}}, \code{\link{filtrationDiag}}
}

\examples{
# input data generated from a circle
X <- circleUnif(n = 10)

# alpha complex filtration
FltAlphaComplex <- alphaComplexFiltration(X = X, printProgress = TRUE)
}

\keyword{methods}
