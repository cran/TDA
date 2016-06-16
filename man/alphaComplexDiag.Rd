\name{alphaComplexDiag}

\alias{alphaComplexDiag}

\title{
  Alpha Complex Persistence Diagram
}

\description{
  The function \code{alphaComplexDiag} computes the persistence diagram of the alpha complex filtration built on top of a point cloud.
}

\usage{
alphaComplexDiag(X, library = "GUDHI", printProgress = FALSE)
}

\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates, used by the function \code{FUN}, where \eqn{n} is the number of points stored in \code{X} and \eqn{d} is the dimension of the space.
  }

  \item{library}{
    Currently only the library "GUDHI" can be used to compute the persistence diagram.
  }

  \item{printProgress}{
    if \code{TRUE}, a progress bar is printed. The default value is \code{FALSE}.
  }
}

\details{
  The function \code{alphaComplexDiag} constructs the alpha complex filtration, using the C++ library \href{ https://project.inria.fr/gudhi/software }{GUDHI}.
  Currently the function \code{alphaComplexDiag} only works on 64-bit R.
  See refereneces.
}

\value{
  The function \code{alphaComplexDiag} returns a list with the following elements:
  \item{diagram}{
    an object of class \code{diagram}, a \eqn{P} by 3 matrix, where \eqn{P} is the number of points in the resulting persistence diagram. The first column stores the dimension of each feature (0 for components, 1 for loops, 2 for voids, etc). Second and third columns are Birth and Death of the features.
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
  \code{\link{summary.diagram}}, \code{\link{plot.diagram}}, \code{\link{alphaShapeDiag}}, \code{\link{gridDiag}}, \code{\link{ripsDiag}}
}

\examples{
# input data generated from a circle
X <- circleUnif(n = 30)

# persistence diagram of alpha complex
DiagAlphaComplex <- alphaComplexDiag(X = X, printProgress = TRUE)

# plot
plot(DiagAlphaComplex[["diagram"]])
}

\keyword{methods}
