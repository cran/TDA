\name{gridFiltration}

\alias{gridFiltration}

\title{
  Persistence Diagram of a function over a Grid
}

\description{
  The function \code{gridFiltration} computes the Persistence Diagram of a filtration of sublevel sets (or superlevel sets) of a function evaluated over a grid of points in arbitrary dimension \code{d}.
}

\usage{
gridFiltration(
    X = NULL, FUN = NULL, lim = NULL, by = NULL, FUNvalues = NULL,
    maxdimension = max(NCOL(X), length(dim(FUNvalues))) - 1,
    sublevel = TRUE, printProgress = FALSE, ...)
}

\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates, used by the function \code{FUN}, where \eqn{n} is the number of points stored in \code{X} and \eqn{d} is the dimension of the space. NULL if this option is not used. The default value is \code{NULL}.
  }

  \item{FUN}{
    a function whose inputs are 1) an \eqn{n} by \eqn{d} matrix of coordinates \code{X}, 2) an \eqn{m} by \eqn{d} matrix of coordinates \code{Grid}, 3) an optional smoothing parameter, and returns a numeric vector of length \eqn{m}. For example see \code{\link{distFct}},  \code{\link{kde}}, and \code{\link{dtm}}  which compute the distance function, the kernel density estimator and the distance to measure, over a grid of points using the input \code{X}. Note that \code{Grid} is not an input of \code{gridFiltration}, but is automatically computed by the function using \code{lim}, and \code{by}. NULL if this option is not used. The default value is \code{NULL}.
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

  \item{printProgress}{
    if \code{TRUE} a progress bar is printed. The default value is \code{FALSE}.
  }

  \item{...}{
    additional parameters for the function \code{FUN}.
  }
}

\details{
  If the values of \code{X}, \code{FUN} are set, then \code{FUNvalues} should be \code{NULL}. In this case, \code{gridFiltration} evaluates the function \code{FUN} over a grid. If the value of \code{FUNvalues} is set, then \code{X}, \code{FUN} should be \code{NULL}. In this case, \code{FUNvalues} is used as function values over the grid.

  Once function values are either computed or given, \code{gridFiltration} constructs a filtration by triangulating the grid and considering the simplices determined by the values of the function of dimension up to \code{maxdimension+1}.
}

\value{
  The function \code{gridFiltration} returns a list with the following elements:

  \item{cmplx}{
    a list representing the complex, same as \code{cmplx} in the input but sorted by the filtration values.
  }

  \item{values}{
    a vector representing the filtration values.
  }

  \item{increasing}{
    a logical variable indicating if the filtration values are in increasing order (\code{TRUE}) or in decreasing order (\code{FALSE}).
  }
}

\references{
  Fasy B, Lecci F, Rinaldo A, Wasserman L, Balakrishnan S, Singh A (2013). "Statistical Inference For Persistent Homology." (arXiv:1303.7117). Annals of Statistics.

  Morozov D (2007). "Dionysus, a C++ library for computing persistent homology." \url{ http://www.mrzv.org/software/dionysus/ }

  Bauer U, Kerber M, Reininghaus J (2012). "PHAT, a software library for persistent homology." \url{ https://bitbucket.org/phat-code/phat/ }
}

\author{
  Brittany T. Fasy, Jisu Kim, and Fabrizio Lecci
}

\note{
  The user can decide to use either the C++ library \href{ https://project.inria.fr/gudhi/software/ }{GUDHI}, \href{ http://www.mrzv.org/software/dionysus/ }{Dionysus}, or \href{ https://bitbucket.org/phat-code/phat/ }{PHAT}. See references.

  Since dimension of simplicial complex from grid points in \eqn{R^d} is up to \eqn{d}, homology of dimension \eqn{\ge d} is trivial. Hence setting \code{maxdimension} with values \eqn{\ge d} is equivalent to \code{maxdimension=d-1}.
}

\seealso{
  \code{\link{summary.diagram}}, \code{\link{plot.diagram}},
  \code{\link{distFct}}, \code{\link{kde}}, \code{\link{kernelDist}}, \code{\link{dtm}},
  \code{\link{alphaComplexDiag}}, \code{\link{alphaComplexDiag}}, \code{\link{ripsDiag}}
}

\examples{
# input data
n <- 10
XX <- circleUnif(n)

## Ranges of the grid
Xlim <- c(-1, 1)
Ylim <- c(-1, 1)
lim <- cbind(Xlim, Ylim)
by <- 1

#Distance Function Diagram of the sublevel sets
FltGrid <- gridFiltration(
  XX, distFct, lim = lim, by = by, sublevel = TRUE, printProgress = TRUE) 

}

\keyword{methods}
