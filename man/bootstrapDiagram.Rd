\name{bootstrapDiagram}

\alias{bootstrapDiagram}

\title{
  Bootstrapped Confidence Set for a Persistence Diagram, using the Bottleneck Distance (or the Wasserstein distance).
}

\description{
  The function \code{bootstrapDiagram} computes a \code{(1-alpha)} confidence set for the Persistence Diagram of a filtration of sublevel sets (or superlevel sets) of a function evaluated over a grid of points. The function returns the (\code{1-alpha}) quantile of \code{B} bottleneck distances (or Wasserstein distances), computed in \code{B} iterations of the bootstrap algorithm.
}

\usage{
bootstrapDiagram(X, FUN, lim, by, maxdimension = length(lim) / 2 - 1,
                 sublevel = TRUE, library = "Dionysus", B = 30, alpha = 0.05,
                 distance = "bottleneck", dimension = min(1, maxdimension),
                 p = 1, parallel = FALSE, printProgress = FALSE, weight = NULL, ...)
}

\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates, used by the function \code{FUN}, where \eqn{n} is the number of points stored in \code{X} and \eqn{d} is the dimension of the space.
  }

  \item{FUN}{
    a function whose inputs are 1) an \eqn{n} by \eqn{d} matrix of coordinates \code{X}, 2) an \eqn{m} by \eqn{d} matrix of coordinates \code{Grid}, 3) an optional smoothing parameter, and returns a numeric vector of length \eqn{m}. For example see \code{\link{distFct}},  \code{\link{kde}}, and \code{\link{dtm}}  which compute the distance function, the kernel density estimator and the distance to measure, over a grid of points using the input \code{X}. Note that \code{Grid} is not an input of \code{bootstrapDiagram}, but is automatically computed by the function using \code{lim} and \code{by}. 
  }  

  \item{lim}{
    a \eqn{2} by \eqn{d} matrix, where each column specifies the range of each dimension of the grid, over which the function \code{FUN} is evaluated.
  }

  \item{by}{
    either a number or a vector of length \eqn{d} specifying space between points of the grid in each dimension. If a number is given, then same space is used in each dimension.
  }

  \item{maxdimension}{
    a number that indicates the maximum dimension to compute persistent homology to. The default value is \eqn{d - 1}, which is (dimension of embedding space - 1).
  }

  \item{sublevel}{
    a logical variable indicating if the Persistence Diagram should be computed for sublevel sets (\code{TRUE}) or superlevel sets (\code{FALSE}) of the function. The default value is \code{TRUE}.
  }

  \item{library}{
    The user can compute the persistence diagram using either the library \code{"Dionysus"}, or \code{"PHAT"}. The default value is \code{"Dionysus"}.
  }

  \item{B}{
    the number of bootstrap iterations. The default value is \code{30}.
  }

  \item{alpha}{
    The function \code{bootstrapDiagram} returns a (\code{1 - alpha}) quantile. The default value is \code{0.05}.
  }

  \item{distance}{
    a string specifying the distance to be used for persistence diagrams: either \code{"bottleneck"} or \code{"wasserstein"}. The default value is \code{"bottleneck"}.
  }

  \item{dimension}{
    \code{dimension} is an integer or a vector specifying the dimension of the features used to compute the bottleneck distance. \code{0} for connected components, \code{1} for loops, \code{2} for voids, and so on. The default value is \code{1} if \eqn{maxdimension \ge 1}, and else \code{0}.
  }

  \item{p}{
    if \code{distance == "wasserstein"}, then \code{p} is an integer specifying the power to be used in the computation of the Wasserstein distance. The default value is \code{1}.
  }

  \item{parallel}{
    logical: if \code{TRUE} the bootstrap iterations are parallelized, using the library \code{parallel}. The default value is \code{FALSE}.
  }

  \item{printProgress}{
    if \code{TRUE} a progress bar is printed. The default value is \code{FALSE}.
  }

  \item{weight}{
    either NULL, a number, or a vector of length \eqn{n}. If it is NULL, weight is not used. If it is a number, then same weight is applied to each points of \code{X}. If it is a vector, \code{weight} represents weights of each points of \code{X}. The default value is \code{NULL}.
  }

  \item{...}{
    additional parameters for the function \code{FUN}.
  }
}

\details{
  The function \code{bootstrapDiagram} uses \code{\link{gridDiag}} to compute the persistence diagram of the input function using the entire sample. Then the bootstrap algorithm, for \code{B} times, computes the bottleneck distance between the original persistence diagram and the one computed using a subsample. Finally the (\code{1-alpha}) quantile of these \code{B} values is returned. See (Chazal, Fasy, Lecci, Michel, Rinaldo, and Wasserman, 2014) for discussion of the method.
}

\value{
  The function \code{bootstrapDiagram} returns the (\code{1-alpha}) quantile of the values computed by the bootstrap algorithm. 
}

\references{
  Chazal F, Fasy BT, Lecci F, Michel B, Rinaldo A, Wasserman L (2014). "Robust Topological Inference: Distance-To-a-Measure and Kernel Distance." Technical Report. 

  Wasserman L (2004), "All of statistics: a concise course in statistical inference." Springer.

  Morozov D (2007). "Dionysus, a C++ library for computing persistent homology." \url{ http://www.mrzv.org/software/dionysus/ }
}

\author{
  Jisu Kim and Fabrizio Lecci
}

\note{
  The function \code{bootstrapDiagram} uses the C++ library \href{ http://www.mrzv.org/software/dionysus }{Dionysus} for the computation of bottleneck and Wasserstein distances. See references.
}

\seealso{
  \code{\link{bottleneck}}, \code{\link{bootstrapBand}},
  \code{\link{distFct}}, \code{\link{kde}}, \code{\link{kernelDist}}, \code{\link{dtm}},
  \code{\link{summary.diagram}}, \code{\link{plot.diagram}}
}

\examples{
## confidence set for the Kernel Density Diagram

# input data
n <- 400
XX <- circleUnif(n)

## Ranges of the grid
Xlim <- c(-1.8, 1.8)
Ylim <- c(-1.6, 1.6)
lim <- cbind(Xlim, Ylim)
by <- 0.05

h <- .3  #bandwidth for the function kde

#Kernel Density Diagram of the superlevel sets
Diag <- gridDiag(XX, kde, lim = lim, by = by, sublevel = FALSE,
                 printProgress = TRUE, h = h) 

# confidence set
B <- 10       ## the number of bootstrap iterations should be higher!
              ## this is just an example
alpha <- 0.05

cc <- bootstrapDiagram(XX, kde, lim = lim, by = by, sublevel = FALSE, B = B,
          alpha = alpha, dimension = 1, printProgress = TRUE, h = h)

plot(Diag[["diagram"]], band = 2 * cc)
}

\keyword{nonparametric}
\keyword{htest}
