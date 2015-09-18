\name{gridDiag}

\alias{gridDiag}

\title{
  Persistence Diagram of a function over a Grid
}

\description{
  \code{gridDiag} computes the Persistence Diagram of a filtration of sublevel sets (or superlevel sets) of a function evaluated over a grid of points in arbitrary dimension \code{d}.
}

\usage{
gridDiag(X = NULL, FUN = NULL, lim = NULL, by = NULL, FUNvalues = NULL,
         maxdimension = max(NCOL(X), length(dim(FUNvalues))) - 1,
         sublevel = TRUE, library = "Dionysus", location = FALSE,
         printProgress = FALSE, diagLimit = NULL, ...)
}

\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates, used by the function \code{FUN}, where \eqn{n} is the number of points stored in \code{X} and \eqn{d} is the dimension of the space. NULL if this option is not used.
  }

  \item{FUN}{
    a function whose inputs are 1) an \eqn{n} by \eqn{d} matrix of coordinates \code{X}, 2) an \eqn{m} by \eqn{d} matrix of coordinates \code{Grid}, 3) an optional smoothing parameter, and returns a numeric vector of length \eqn{m}. For example see \code{\link{distFct}},  \code{\link{kde}}, and \code{\link{dtm}}  which compute the distance function, the kernel density estimator and the distance to measure, over a grid of points using the input \code{X}. Note that \code{Grid} is not an input of \code{gridDiag}, but is automatically computed by the function using \code{lim}, and \code{by}. NULL if this option is not used.
  }

  \item{lim}{
    a \eqn{2} by \eqn{d} matrix, where each column specifying the range of each dimension of the grid, over which the function \code{FUN} is evaluated. NULL if this option is not used.
  }

  \item{by}{
    either a number or a vector of length \eqn{d} specifying space between points of the grid in each dimension. If a number is given, then same space is used in each dimension. NULL if this option is not used.
  }

  \item{FUNvalues}{
    an \eqn{m1 * m2 * ... * md} array of function values over \eqn{m1 * m2 * ... * md} grid, where \eqn{mi} is the number of scales of grid on \eqn{ith} dimension. NULL if this option is not used.
  }

  \item{maxdimension}{
    a number that indicates the maximum dimension of the homological features to compute: 0 for connected components, 1 for loops, 2 for voids and so on. Default is \eqn{d-1}, which is (dimension of embedding space - 1).
  }

  \item{sublevel}{
    a logical variable indicating if the Persistence Diagram should be computed for sublevel sets (\code{TRUE}) or superlevel sets (\code{FALSE}) of the function. Default is \code{TRUE}.
  }

  \item{library}{
    The user can compute the persistence diagram using either the library 'Dionysus', or 'PHAT'. Default is 'Dionysus'.
  }

  \item{location}{
    if \code{TRUE}, location of birth point and death point of each homological feature is returned. Additionaly if \code{library="Dionysus"}, location of representative cycles of each homological feature is also returned.
  }

  \item{printProgress}{
    if \code{TRUE} a progress bar is printed. Default is \code{FALSE}.
  }

  \item{diagLimit}{
    a number that replaces \code{Inf} (if \code{sublevel} is \code{TRUE}) or \code{-Inf} (if \code{sublevel} is \code{FALSE}) in the Death value of the most persistent connected component. Deafult is \code{NULL} and the max/min of the function is used.
  }

  \item{...}{
    additional parameters for the function \code{FUN}.
  }
}

\details{
  If the values of \code{X}, \code{FUN}, \code{lim}, and \code{by} are set, then \code{FUNvalues} should be \code{NULL}. In this case, \code{gridDiag} evaluates the function \code{FUN} over a grid. If the value of \code{FUNvalues} is set, then \code{X}, \code{FUN}, \code{lim}, and \code{by} should be \code{NULL}. In this case, \code{FUNvalues} is used as function values over the grid.

  Once function values are either computed or given, \code{gridDiag} constructs a filtration by triangulating the grid and considering the simplices determined by the values of the function of dimension up to \code{maxdimension+1}.
}

\value{
  \code{gridDiag} returns a list with the following components:

  \item{diagram}{
    an object of class \code{diagram}, a \eqn{P} by 3 matrix, where \eqn{P} is the number of points in the resulting persistence diagram. The first column stores the dimension of each feature (0 for components, 1 for loops, 2 for voids, etc). Second and third columns are Birth and Death of the features, in case of a filtration constructed using sublevel sets (from -Inf to Inf), or Death and Birth of features, in case of a filtration constructed using superlevel sets (from Inf to -Inf).
  }
  \item{birthLocation}{
    only if \code{location = TRUE}: a \eqn{P} by \eqn{d} matrix, where \eqn{P} is the number of points in the resulting persistence diagram. Each row represents the location of the grid point completing the simplex that gives birth to an homological feature.
  }
  \item{deathLocation}{
    only if \code{location = TRUE}: a \eqn{P} by \eqn{d} matrix, where \eqn{P} is the number of points in the resulting persistence diagram. Each row represents the location of the grid point completing the simplex that kills an homological feature.
  }
  \item{cycleLocation}{
    only if \code{location = TRUE} and \code{library = "Dionysus"}: a list of length \eqn{P}, where \eqn{P} is the number of points in the resulting persistence diagram. Each element is a \eqn{P_i} by \eqn{d} matrix and represents location of \eqn{P_i} grid points on a representative cycle of each homological feature.
  }
}

\references{
  Brittany Fasy, Fabrizio Lecci, Alessandro Rinaldo, Larry Wasserman, Sivaraman Balakrishnan, and Aarti Singh. (2013), "Statistical Inference For Persistent Homology", (arXiv:1303.7117). To appear, Annals of Statistics.

  Dmitriy Morozov, "Dionysus, a C++ library for computing persistent homology". http://www.mrzv.org/software/dionysus/

  Ulrich Bauer, Michael Kerber, Jan Reininghaus, "PHAT, a software library for persistent homology". https://code.google.com/p/phat/
}

\author{
  Brittany T. Fasy, Jisu Kim, and Fabrizio Lecci
}

\note{
  The user can decide to use either the C++ library Dionysus or the C++ library PHAT. See references.

  Since dimension of simplicial complex from grid points in \eqn{R^d} is up to \eqn{d}, homology of dimension \eqn{>=d} is trivial. Hence setting \code{maxdimension} with values \eqn{>=d} is equivalent to \code{maxdimension=d-1}.
}

\seealso{
  \code{\link{summary.diagram}}, \code{\link{plot.diagram}}, \code{\link{distFct}}, \code{\link{kde}}, \code{\link{kernelDist}}, \code{\link{dtm}}, \code{\link{ripsDiag}}
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
                  location = TRUE, printProgress = TRUE, h = h) 
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
for (i in seq(along = one))
{
  points(Diag2[["birthLocation"]][one[i], ], pch = 15, cex = 3, col = i)
  points(Diag2[["deathLocation"]][one[i], ], pch = 17, cex = 3, col = i)
  points(Diag2[["cycleLocation"]][[one[i]]], pch = 19, cex = 1, col = i)
}
}

\keyword{methods}
