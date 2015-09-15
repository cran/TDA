\name{wasserstein}

\alias{wasserstein}

\title{
  Wasserstein distance between two persistence diagrams
}

\description{
  This function computes the Wasserstein distance between two persistence diagrams.
}

\usage{
  wasserstein(Diag1, Diag2, p = 1, dimension = 1)
}

\arguments{
  \item{Diag1}{
    an object of class \code{diagram} or a matrix (\eqn{n} by 3) that stores dimension, birth and death of \eqn{n} topological features.
  }

  \item{Diag2}{
    an object of class \code{diagram} or a matrix (\eqn{m} by 3) that stores dimension, birth and death of \eqn{m} topological features.
  }

  \item{p}{
    integer specifying the power to be used in the computation of the Wasserstein distance. Default is 1.
  }

  \item{dimension}{
    an integer or a vector specifying the dimension of the features used to compute the wasserstein distance. 0 for connected components, 1 for loops, 2 for voids and so on. Default is 1 (loops).
  }
}

\details{
  The Wasserstein distance between two diagrams is the cost of the optimal matching between points of the two diagrams. When a vector is given for \code{dimension}, then maximum among bottleneck distances using each element in \code{dimension} is returned. This function is an R wrapper of the function "wasserstein_distance" in the C++ library Dionysus. See references.
}

\value{
  Returns the value of the Wasserstein distance between the two persistence diagrams.
}

\references{
  Dmitriy Morozov, "Dionysus, a C++ library for computing persistent homology". http://www.mrzv.org/software/dionysus/

  Herbert Edelsbrunner and John Harer (2010), Computational topology: an introduction. American Mathematical Society.
}

\author{
Jisu Kim, Fabrizio Lecci
}

\seealso{
  \code{\link{bottleneck}}, \code{\link{ripsDiag}},  \code{\link{gridDiag}},
  \code{\link{plot.diagram}}
}

\examples{
XX1 <- circleUnif(20)
XX2 <- circleUnif(20, r = 0.2)

DiagLim <- 5
maxdimension <- 1

Diag1 <- ripsDiag(XX1, maxdimension, DiagLim, printProgress = FALSE)
Diag2 <- ripsDiag(XX2, maxdimension, DiagLim, printProgress = FALSE)

wassersteinDist <- wasserstein(Diag1[["diagram"]], Diag2[["diagram"]], p = 1,
                               dimension = 1)
print(wassersteinDist)
}

\keyword{methods}
\keyword{optimize}
