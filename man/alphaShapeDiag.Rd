\name{alphaShapeDiag}

\alias{alphaShapeDiag}

\title{
  Persistence Diagram of Alpha Shape in 3d
}

\description{
  The function \code{alphaShapeDiag} computes the persistence diagram of the alpha shape filtration built on top of a point cloud in 3 dimension.
}

\usage{
alphaShapeDiag(
    X, maxdimension = NCOL(X) - 1, library = "GUDHI", location = FALSE,
    printProgress = FALSE)
}

\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates, used by the function \code{FUN}, where \eqn{n} is the number of points stored in \code{X} and \eqn{d} is the dimension of the space.
    Currently \eqn{d} should be 3.
  }

  \item{maxdimension}{
    integer: max dimension of the homological features to be computed. (e.g. 0 for connected components, 1 for connected components and loops, 2 for connected components, loops, voids, etc.)
  }

  \item{library}{
    either a string or a vector of length two. When a vector is given, the first element specifies which library to compute the Alpha Shape filtration, and the second element specifies which library to compute the persistence diagram. If a string is used, then the same library is used. For computing the Alpha Shape filtration, the user can use the library \code{"GUDHI"}, and is also the default value. For computing the persistence diagram, the user can choose either the library \code{"GUDHI"}, \code{"Dionysus"}, or \code{"PHAT"}. The default value is \code{"GUDHI"}.
  }

  \item{location}{
    if \code{TRUE} and if \code{"Dionysus"} or \code{"PHAT"} is used for computing the persistence diagram, location of birth point and death point of each homological feature is returned. Additionaly if \code{library="Dionysus"}, location of representative cycles of each homological feature is also returned. The default value is \code{FALSE}.
  }

  \item{printProgress}{
    if \code{TRUE}, a progress bar is printed. The default value is \code{FALSE}.
  }
}

\details{
  The function \code{alphaShapeDiag} constructs the Alpha Shape filtration, using the C++ library \href{ https://project.inria.fr/gudhi/software }{GUDHI}.
  Then for computing the persistence diagram from the Alpha Shape filtration, the user can use either the C++ library \href{ https://project.inria.fr/gudhi/software/ }{GUDHI}, \href{ http://www.mrzv.org/software/dionysus/ }{Dionysus}, or \href{ https://bitbucket.org/phat-code/phat/ }{PHAT}.
  See refereneces.
}

\value{
  The function \code{alphaShapeDiag} returns a list with the following elements:
  \item{diagram}{
    an object of class \code{diagram}, a \eqn{P} by 3 matrix, where \eqn{P} is the number of points in the resulting persistence diagram. The first column stores the dimension of each feature (0 for components, 1 for loops, 2 for voids, etc). Second and third columns are Birth and Death of the features.
  }
  \item{birthLocation}{
    only if \code{location=TRUE} and if \code{"Dionysus"} or \code{"PHAT"} is used for computing the persistence diagram: a \eqn{P} by \eqn{d} matrix, where \eqn{P} is the number of points in the resulting persistence diagram. Each row represents the location of the grid point completing the simplex that gives birth to an homological feature.
  }
  \item{deathLocation}{
    only if \code{location=TRUE} and if \code{"Dionysus"} or \code{"PHAT"} is used for computing the persistence diagram: a \eqn{P} by \eqn{d} matrix, where \eqn{P} is the number of points in the resulting persistence diagram. Each row represents the location of the grid point completing the simplex that kills an homological feature.
  }
  \item{cycleLocation}{
    only if \code{location=TRUE} and if \code{"Dionysus"} is used for computing the persistence diagram: a list of length \eqn{P}, where \eqn{P} is the number of points in the resulting persistence diagram. Each element is a \eqn{P_i} by \eqn{h_i +1} by \eqn{d} array for \eqn{h_i} dimensional homological feature. It represents location of \eqn{h_i +1} vertices of \eqn{P_i} simplices, where \eqn{P_i} simplices constitutes the \eqn{h_i} dimensional homological feature.
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
X <- sphereUnif(n = 500, d = 2)

# persistence diagram of alpha shape
DiagAlphaShape <-
    alphaShapeDiag(X = X, library = "GUDHI", printProgress = TRUE)

# plot
plot(DiagAlphaShape[["diagram"]])

}

\keyword{methods}
