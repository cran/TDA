\name{filtrationDiag}

\alias{filtrationDiag}

\title{
  Persistence Diagram of Filtration
}

\description{
  The function \code{filtrationDiag} computes the persistence diagram of the filtration.
}

\usage{
filtrationDiag(
    filtration, maxdimension, library = "GUDHI", location = FALSE,
    printProgress = FALSE, diagLimit = NULL)
}

\arguments{
  \item{filtration}{
    a list representing the input filtration. This list consists of three components: \code{"cmplx"}, a list representing the complex, \code{"values"}, a vector representing the filtration values, and \code{"increasing"}, a logical variable indicating if the filtration values are in increasing order or in decreasing order.
  }

  \item{maxdimension}{
    integer: max dimension of the homological features to be computed. (e.g. 0 for connected components, 1 for connected components and loops, 2 for connected components, loops, voids, etc.)
  }

  \item{library}{
    a string specifying which library to compute the persistence diagram. The user can choose either the library \code{"GUDHI"} or \code{"Dionysus"}. The default value is \code{"GUDHI"}.
  }

  \item{location}{
    if \code{TRUE} and if \code{"Dionysus"} is used for computing the persistence diagram, location of birth point, death point, and representative cycles, of each homological feature is returned.
  }

  \item{printProgress}{
    logical: if \code{TRUE}, a progress bar is printed. The default value is \code{FALSE}.
  }

  \item{diagLimit}{
    a number that replaces \code{Inf} in the persistence diagram. The default value is \code{NULL} and \code{Inf} value in the persistence diagram will not be replaced.
  }
}

\details{
  The user can decide to use either the C++ library \href{ https://project.inria.fr/gudhi/software }{GUDHI} or \href{ http://www.mrzv.org/software/dionysus }{Dionysus}.
  See refereneces.
}

\value{
  The function \code{filtrationDiag} returns a list with the following elements:

  \item{diagram}{
    an object of class \code{diagram}, a \eqn{P} by 3 matrix, where \eqn{P} is the number of points in the resulting persistence diagram. The first column contains the dimension of each feature (0 for components, 1 for loops, 2 for voids, etc.). Second and third columns are Birth and Death of the features.
  }

  \item{birthLocation}{
    only if \code{location=TRUE} and if \code{"Dionysus"} is used for computing the persistence diagram: a vector of length \eqn{P}. Each row represents the index of the vertex completing the simplex that gives birth to an homological feature.
  }

  \item{deathLocation}{
    only if \code{location=TRUE} and if \code{"Dionysus"} is used for computing the persistence diagram: a vector of length \eqn{P}. Each row represents the index of the vertex completing the simplex that kills an homological feature.
  }

  \item{cycleLocation}{
    only if \code{location=TRUE} and if \code{"Dionysus"} is used for computing the persistence diagram: a \eqn{P_i} by \eqn{h_i +1} matrix for \eqn{h_i} dimensional homological feature. It represents index of \eqn{h_i +1} vertices of \eqn{P_i} simplices on a representative cycle of the \eqn{h_i} dimensional homological feature.
  }
}

\references{
  Maria C (2014). "GUDHI, Simplicial Complexes and Persistent Homology Packages." \url{ https://project.inria.fr/gudhi/software/ }.

  Morozov D (2007). "Dionysus, a C++ library for computing persistent homology". \url{ http://www.mrzv.org/software/dionysus/ }

  Edelsbrunner H, Harer J (2010). "Computational topology: an introduction." American Mathematical Society.

  Fasy B, Lecci F, Rinaldo A, Wasserman L, Balakrishnan S, Singh A (2013). "Statistical Inference For Persistent Homology." (arXiv:1303.7117). Annals of Statistics.
}

\author{
  Jisu Kim
}

\seealso{
  \code{\link{summary.diagram}}, \code{\link{plot.diagram}}
}

\examples{
n <- 5
X <- cbind(cos(2*pi*seq_len(n)/n), sin(2*pi*seq_len(n)/n))
maxdimension <- 1
maxscale <- 1.5
dist <- "euclidean"
library <- "Dionysus"

FltRips <- ripsFiltration(X = X, maxdimension = maxdimension,
               maxscale = maxscale, dist = "euclidean", library = "Dionysus",
               printProgress = TRUE)

DiagFltRips <- filtrationDiag(filtration = FltRips, maxdimension = maxdimension,
                   library = "Dionysus", location = TRUE, printProgress = TRUE)

plot(DiagFltRips[["diagram"]])


FUNvalues <- X[, 1] + X[, 2]

FltFun <- funFiltration(FUNvalues = FUNvalues, cmplx = FltRips[["cmplx"]])

DiagFltFun <- filtrationDiag(filtration = FltFun, maxdimension = maxdimension,
                             library = "Dionysus", location = TRUE, printProgress = TRUE)

plot(DiagFltFun[["diagram"]], diagLim = c(-2, 5))
}
\keyword{methods}
