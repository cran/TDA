\name{ripsDiag}

\alias{ripsDiag}

\title{
  Rips Persistence Diagram
}

\description{
  The function \code{ripsDiag} computes the persistence diagram of the Rips filtration built on top of a point cloud.
}

\usage{
ripsDiag(
    X, maxdimension, maxscale, dist = "euclidean",
    library = "GUDHI", location = FALSE, printProgress = FALSE)
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
    either a string or a vector of length two. When a vector is given, the first element specifies which library to compute the Rips filtration, and the second element specifies which library to compute the persistence diagram. If a string is used, then the same library is used. For computing the Rips filtration, if \code{dist = "euclidean"}, the user can use either the library \code{"GUDHI"} or \code{"Dionysus"}. If \code{dist = "arbitrary"}, the user can use either the library \code{"Dionysus"}. The default value is \code{"GUDHI"} if \code{dist = "euclidean"}, and \code{"Dionysus"} if \code{dist == "arbitrary"}. When \code{"GUDHI"} is used for \code{dist = "arbitrary"}, \code{"Dionysus"} is implicitly used. For computing the persistence diagram, the user can choose either the library \code{"GUDHI"}, \code{"Dionysus"}, or \code{"PHAT"}. The default value is \code{"GUDHI"}.
  }

  \item{location}{
    if \code{TRUE} and if \code{"Dionysus"} or \code{"PHAT"} is used for computing the persistence diagram, location of birth point and death point of each homological feature is returned. Additionaly if \code{library="Dionysus"}, location of representative cycles of each homological feature is also returned.
  }

  \item{printProgress}{
    logical: if \code{TRUE}, a progress bar is printed. The default value is \code{FALSE}.
  }
}

\details{
  For Rips filtration based on Euclidean distance of the input point cloud, the user can decide to use either the C++ library \href{ https://project.inria.fr/gudhi/software }{GUDHI} or \href{ http://www.mrzv.org/software/dionysus }{Dionysus}.
  For Rips filtration based on arbitrary distance, the user can decide to the C++ library \href{ http://www.mrzv.org/software/dionysus }{Dionysus}.
  Then for computing the persistence diagram from the Rips filtration, the user can use either the C++ library \href{ https://project.inria.fr/gudhi/software/ }{GUDHI}, \href{ http://www.mrzv.org/software/dionysus/ }{Dionysus}, or \href{ https://bitbucket.org/phat-code/phat/ }{PHAT}.
  See refereneces.
}

\value{
  The function \code{ripsDiag} returns a list with the following elements:

  \item{diagram}{
    an object of class \code{diagram}, a \eqn{P} by 3 matrix, where \eqn{P} is the number of points in the resulting persistence diagram. The first column contains the dimension of each feature (0 for components, 1 for loops, 2 for voids, etc.). Second and third columns are Birth and Death of the features.
  }

  \item{birthLocation}{
    only if \code{location=TRUE} and if \code{"Dionysus"} or \code{"PHAT"} is used for computing the persistence diagram: if \code{dist="euclidean"}, then \code{birthLocation} is a \eqn{P} by \eqn{d} matrix, where \eqn{P} is the number of points in the resulting persistence diagram. Each row represents the location of the data point completing the simplex that gives birth to an homological feature. If \code{dist="arbitrary"}, then \code{birthLocation} is a vector of length \eqn{P}. Each row represents the index of the data point completing the simplex that gives birth to an homological feature.
  }

  \item{deathLocation}{
    only if \code{location=TRUE} and if \code{"Dionysus"} or \code{"PHAT"} is used for computing the persistence diagram: if \code{dist="euclidean"}, then \code{deathLocation} is a \eqn{P} by \eqn{d} matrix, where \eqn{P} is the number of points in the resulting persistence diagram. Each row represents the location of the data point completing the simplex that kills an homological feature. If \code{dist="arbitrary"}, then \code{deathLocation} is a vector of length \eqn{P}. Each row represents the index of the data point completing the simplex that kills an homological feature.
  }

  \item{cycleLocation}{
    only if \code{location=TRUE} and if \code{"Dionysus"} is used for computing the persistence diagram: if \code{dist="euclidean"}, then \code{cycleLocation} is a list of length \eqn{P}, where \eqn{P} is the number of points in the resulting persistence diagram. Each element is a \eqn{P_i} by \eqn{h_i +1} by \eqn{d} array for \eqn{h_i} dimensional homological feature. It represents location of \eqn{h_i +1} vertices of \eqn{P_i} simplices, where \eqn{P_i} simplices constitutes the \eqn{h_i} dimensional homological feature. If \code{dist = "arbitrary"}, then each element is a \eqn{P_i} by \eqn{h_i +1} matrix for  for \eqn{h_i} dimensional homological feature. It represents index of \eqn{h_i +1} vertices of \eqn{P_i} simplices on a representative cycle of the \eqn{h_i} dimensional homological feature.
  }
}

\references{
  Maria C (2014). "GUDHI, Simplicial Complexes and Persistent Homology Packages." \url{ https://project.inria.fr/gudhi/software/ }.

  Morozov D (2007). "Dionysus, a C++ library for computing persistent homology". \url{ http://www.mrzv.org/software/dionysus/ }

  Edelsbrunner H, Harer J (2010). "Computational topology: an introduction." American Mathematical Society.

  Fasy B, Lecci F, Rinaldo A, Wasserman L, Balakrishnan S, Singh A (2013). "Statistical Inference For Persistent Homology." (arXiv:1303.7117). Annals of Statistics.
}

\author{
  Brittany T. Fasy, Jisu Kim, Fabrizio Lecci, and Clement Maria
}

\seealso{
  \code{\link{summary.diagram}}, \code{\link{plot.diagram}}, \code{\link{gridDiag}}
}

\examples{
## EXAMPLE 1: rips diagram for circles (euclidean distance)
XX <- circleUnif(30)
maxscale <- 5
maxdimension <- 1
## note that the input XX is a point cloud
Diag <- ripsDiag(XX, maxdimension, maxscale, library = "Dionysus",
    location = TRUE, printProgress = TRUE)

# plot
layout(matrix(c(1, 3, 2, 2), 2, 2))
plot(XX, cex = 0.5, pch = 19)
title(main = "Data")
plot(Diag[["diagram"]])
title(main = "rips Diagram")
one <- which(Diag[["diagram"]][, 1] == 1)
plot(XX, col = 2, main = "Representative loop of data points")
for (i in seq(along = one)) {
  for (j in seq_len(dim(Diag[["cycleLocation"]][[one[i]]])[1])) {
    lines(Diag[["cycleLocation"]][[one[i]]][j, , ], pch = 19, cex = 1, col = i)
  }
}


## EXAMPLE 2: rips diagram with arbitrary distance
## distance matrix for triangle with edges of length: 1,2,4
distX <- matrix(c(0, 1, 2, 1, 0, 4, 2, 4, 0), ncol = 3)
maxscale <- 5
maxdimension <- 1
## note that the input distXX is a distance matrix
DiagTri <- ripsDiag(distX, maxdimension, maxscale, dist = "arbitrary",
                    printProgress = TRUE)
#points with lifetime = 0 are not shown. e.g. the loop of the triangle.
print(DiagTri[["diagram"]])
}
\keyword{methods}
