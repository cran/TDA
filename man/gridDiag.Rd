\name{gridDiag}

\alias{gridDiag}

\title{
  Persistence Diagram of a function over a Grid
}

\description{
  The function \code{gridDiag} computes the Persistence Diagram of a filtration of sublevel sets (or superlevel sets) of a function evaluated over a grid of points in arbitrary dimension \code{d}.
}

\usage{
gridDiag(
    X = NULL, FUN = NULL, lim = NULL, by = NULL, FUNvalues = NULL,
    maxdimension = max(NCOL(X), length(dim(FUNvalues))) - 1,
    sublevel = TRUE, library = "GUDHI", location = FALSE,
    printProgress = FALSE, diagLimit = NULL, ...)
}

\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates, used by the function \code{FUN}, where \eqn{n} is the number of points stored in \code{X} and \eqn{d} is the dimension of the space. NULL if this option is not used. The default value is \code{NULL}.
  }

  \item{FUN}{
    a function whose inputs are 1) an \eqn{n} by \eqn{d} matrix of coordinates \code{X}, 2) an \eqn{m} by \eqn{d} matrix of coordinates \code{Grid}, 3) an optional smoothing parameter, and returns a numeric vector of length \eqn{m}. For example see \code{\link{distFct}},  \code{\link{kde}}, and \code{\link{dtm}}  which compute the distance function, the kernel density estimator and the distance to measure, over a grid of points using the input \code{X}. Note that \code{Grid} is not an input of \code{gridDiag}, but is automatically computed by the function using \code{lim}, and \code{by}. NULL if this option is not used. The default value is \code{NULL}.
  }

  \item{lim}{
    a \eqn{2} by \eqn{d} matrix, where each column specifying the range of each dimension of the grid, over which the function \code{FUN} is evaluated. NULL if this option is not used. The default value is \code{NULL}.
  }

  \item{by}{
    either a number or a vector of length \eqn{d} specifying space between points of the grid in each dimension. If a number is given, then same space is used in each dimension. NULL if this option is not used. The default value is \code{NULL}.
  }

  \item{FUNvalues}{
    an \eqn{m1 * m2 * ... * md} array of function values over \eqn{m1 * m2 * ... * md} grid, where \eqn{mi} is the number of scales of grid on \eqn{ith} dimension. NULL if this option is not used. The default value is \code{NULL}.
  }

  \item{maxdimension}{
    a number that indicates the maximum dimension of the homological features to compute: \code{0} for connected components, \code{1} for loops, \code{2} for voids and so on. The default value is \eqn{d - 1}, which is (dimension of embedding space - 1).
  }

  \item{sublevel}{
    a logical variable indicating if the Persistence Diagram should be computed for sublevel sets (\code{TRUE}) or superlevel sets (\code{FALSE}) of the function. The default value is \code{TRUE}.
  }

  \item{library}{
    a string specifying which library to compute the persistence diagram. The user can choose either the library \code{"GUDHI"}, \code{"Dionysus"}, or \code{"PHAT"}. The default value is \code{"GUDHI"}.
  }

  \item{location}{
    if \code{TRUE} and if \code{"Dionysus"} or \code{"PHAT"} is used for computing the persistence diagram, location of birth point and death point of each homological feature is returned. Additionaly if \code{library="Dionysus"}, location of representative cycles of each homological feature is also returned. The default value is \code{FALSE}.
  }

  \item{printProgress}{
    if \code{TRUE} a progress bar is printed. The default value is \code{FALSE}.
  }

  \item{diagLimit}{
    a number that replaces \code{Inf} (if \code{sublevel} is \code{TRUE}) or \code{-Inf} (if \code{sublevel} is \code{FALSE}) in the Death value of the most persistent connected component. The default value is \code{NULL} and the max/min of the function is used.
  }

  \item{...}{
    additional parameters for the function \code{FUN}.
  }
}

\details{
  If the values of \code{X}, \code{FUN} are set, then \code{FUNvalues} should be \code{NULL}. In this case, \code{gridDiag} evaluates the function \code{FUN} over a grid. If the value of \code{FUNvalues} is set, then \code{X}, \code{FUN} should be \code{NULL}. In this case, \code{FUNvalues} is used as function values over the grid. If \code{location=TRUE}, then \code{lim}, and \code{by} should be set.

  Once function values are either computed or given, \code{gridDiag} constructs a filtration by triangulating the grid and considering the simplices determined by the values of the function of dimension up to \code{maxdimension+1}.
}

\value{
  The function \code{gridDiag} returns a list with the following components:

  \item{diagram}{
    an object of class \code{diagram}, a \eqn{P} by 3 matrix, where \eqn{P} is the number of points in the resulting persistence diagram. The first column stores the dimension of each feature (0 for components, 1 for loops, 2 for voids, etc). Second and third columns are Birth and Death of the features, in case of a filtration constructed using sublevel sets (from -Inf to Inf), or Death and Birth of features, in case of a filtration constructed using superlevel sets (from Inf to -Inf).
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
  Fasy B, Lecci F, Rinaldo A, Wasserman L, Balakrishnan S, Singh A (2013). "Statistical Inference For Persistent Homology." (arXiv:1303.7117). Annals of Statistics.

  Morozov D (2007). "Dionysus, a C++ library for computing persistent homology." \url{ https://www.mrzv.org/software/dionysus/ }

  Bauer U, Kerber M, Reininghaus J (2012). "PHAT, a software library for persistent homology." \url{ https://bitbucket.org/phat-code/phat/ }
}

\author{
  Brittany T. Fasy, Jisu Kim, and Fabrizio Lecci
}

\note{
  The user can decide to use either the C++ library \href{ https://project.inria.fr/gudhi/software/ }{GUDHI}, \href{ https://www.mrzv.org/software/dionysus/ }{Dionysus}, or \href{ https://bitbucket.org/phat-code/phat/ }{PHAT}. See references.

  Since dimension of simplicial complex from grid points in \eqn{R^d} is up to \eqn{d}, homology of dimension \eqn{\ge d} is trivial. Hence setting \code{maxdimension} with values \eqn{\ge d} is equivalent to \code{maxdimension=d-1}.
}

\seealso{
  \code{\link{summary.diagram}}, \code{\link{plot.diagram}},
  \code{\link{distFct}}, \code{\link{kde}}, \code{\link{kernelDist}}, \code{\link{dtm}},
  \code{\link{alphaComplexDiag}}, \code{\link{alphaComplexDiag}}, \code{\link{ripsDiag}}
}

\examples{
## Distance Function Diagram and Kernel Density Diagram

# input data
n <- 300
XX <- circleUnif(n)

## Ranges of the grid
Xlim <- c(-1.8, 1.8)
Ylim <- c(-1.6, 1.6)
lim <- cbind(Xlim, Ylim)
by <- 0.05

h <- .3  #bandwidth for the function kde

#Distance Function Diagram of the sublevel sets
Diag1 <- gridDiag(XX, distFct, lim = lim, by = by, sublevel = TRUE,
                  printProgress = TRUE) 

#Kernel Density Diagram of the superlevel sets
Diag2 <- gridDiag(XX, kde, lim = lim, by = by, sublevel = FALSE,
    library = "Dionysus", location = TRUE, printProgress = TRUE, h = h)
#plot
par(mfrow = c(2, 2))
plot(XX, cex = 0.5, pch = 19)
title(main = "Data")
plot(Diag1[["diagram"]])
title(main = "Distance Function Diagram")
plot(Diag2[["diagram"]])
title(main = "Density Persistence Diagram")
one <- which(Diag2[["diagram"]][, 1] == 1)
plot(XX, col = 2, main = "Representative loop of grid points")
for (i in seq(along = one)) {
  points(Diag2[["birthLocation"]][one[i], , drop = FALSE], pch = 15, cex = 3,
      col = i)
  points(Diag2[["deathLocation"]][one[i], , drop = FALSE], pch = 17, cex = 3,
      col = i)
  for (j in seq_len(dim(Diag2[["cycleLocation"]][[one[i]]])[1])) {
    lines(Diag2[["cycleLocation"]][[one[i]]][j, , ], pch = 19, cex = 1, col = i)
  }
}
}

\keyword{methods}
