\name{ripsFiltration}

\alias{ripsFiltration}

\title{
  Rips Filtration
}

\description{
  The function \code{ripsFiltration} computes the Rips filtration built on top of a point cloud.
}

\usage{
ripsFiltration(
    X, maxdimension, maxscale, dist = "euclidean",
    library = "GUDHI", printProgress = FALSE)
}

\arguments{
  \item{X}{
    If \code{dist="euclidean"}, \code{X} is an \eqn{n} by \eqn{d} matrix of coordinates, where \eqn{n} is the number of points in the \eqn{d}-dimensional euclidean space.
    If \code{dist="arbitrary"}, \code{X} is an \eqn{n} by \eqn{n} matrix of distances of \eqn{n} points.
  }

  \item{maxdimension}{
    integer: max dimension of the homological features to be computed. (e.g. 0 for connected components, 1 for connected components and loops, 2 for connected components, loops, voids, etc.)
  }

  \item{maxscale}{
    number: maximum value of the rips filtration.
  }

  \item{dist}{
    \code{"euclidean"} for Euclidean distance, \code{"arbitrary"} for an arbitrary distance given in input as a distance matrix.
  }

  \item{library}{
    a string specifying which library to compute the Rips filtration. If \code{dist = "euclidean"}, the user can use either the library \code{"GUDHI"} or \code{"Dionysus"}. If \code{dist = "arbitrary"}, the user can use the library \code{"Dionysus"}. The default value is \code{"GUDHI"} if \code{dist = "euclidean"}, and \code{"Dionysus"} if \code{dist == "arbitrary"}. When \code{"GUDHI"} is used for \code{dist = "arbitrary"}, \code{"Dionysus"} is implicitly used.
  }

  \item{printProgress}{
    logical: if \code{TRUE}, a progress bar is printed. The default value is \code{FALSE}.
  }
}

\details{
  For Rips filtration based on Euclidean distance of the input point cloud, the user can decide to use either the C++ library \href{ https://project.inria.fr/gudhi/software }{GUDHI} or \href{ http://www.mrzv.org/software/dionysus }{Dionysus}.
  For Rips filtration based on arbitrary distance, the user can use the C++ library \href{ http://www.mrzv.org/software/dionysus }{Dionysus}.
  See refereneces.
}

\value{
  The function \code{ripsFiltration} returns a list with the following elements:

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
    only if \code{dist = "euclidean"}: a matrix representing the coordinates of vertices. Its i-th row represents the coordinate of i-th vertex.
  }
}

\references{
  Maria C (2014). "GUDHI, Simplicial Complexes and Persistent Homology Packages." \url{ https://project.inria.fr/gudhi/software/ }.

  Morozov D (2007). "Dionysus, a C++ library for computing persistent homology". \url{ http://www.mrzv.org/software/dionysus/ }

  Edelsbrunner H, Harer J (2010). "Computational topology: an introduction." American Mathematical Society.
}

\author{
  Jisu Kim
}

\seealso{
  \code{\link{ripsDiag}}, \code{\link{filtrationDiag}}
}

\examples{
n <- 5
X <- cbind(cos(2*pi*seq_len(n)/n), sin(2*pi*seq_len(n)/n))
maxdimension <- 1
maxscale <- 1.5

FltRips <- ripsFiltration(X = X, maxdimension = maxdimension,
               maxscale = maxscale, dist = "euclidean", library = "GUDHI",
               printProgress = TRUE)
}
\keyword{methods}
