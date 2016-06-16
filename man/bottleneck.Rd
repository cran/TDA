\name{bottleneck}

\alias{bottleneck}

\title{
  Bottleneck distance between two persistence diagrams
}

\description{
  The function \code{bottleneck} computes the bottleneck distance between two persistence diagrams.
}

\usage{
bottleneck(Diag1, Diag2, dimension = 1)
}

\arguments{
  \item{Diag1}{
    an object of class \code{diagram} or a matrix (\eqn{n} by 3) that stores dimension, birth and death of \eqn{n} topological features.
  }

  \item{Diag2}{
    an object of class \code{diagram} or a matrix (\eqn{m} by 3) that stores dimension, birth and death of \eqn{m} topological features.
  }

  \item{dimension}{
    an integer or a vector specifying the dimension of the features used to compute the bottleneck distance. \code{0} for connected components, \code{1} for loops, \code{2} for voids and so on. The default value is \code{1} (loops). The default value is \code{1}.
  }
}

\details{
  The bottleneck distance between two diagrams is the cost of the optimal matching between points of the two diagrams. Note that all the diagonal points are included in the persistence diagrams when computing the optimal matching. When a vector is given for \code{dimension}, then maximum among bottleneck distances using each element in \code{dimension} is returned. The function \code{bottleneck} is an R wrapper of the function "bottleneck_distance" in the C++ library \href{ http://www.mrzv.org/software/dionysus }{Dionysus}. See references.
}

\value{
  The function \code{bottleneck} returns the value of the bottleneck distance between the two persistence diagrams.
}

\references{
  Morozov D (2007). "Dionysus, a C++ library for computing persistent homology." \url{ http://www.mrzv.org/software/dionysus/ }

  Edelsbrunner H, Harer J (2010). "Computational topology: an introduction." American Mathematical Society.
}

\author{
  Jisu Kim and Fabrizio Lecci
}

\seealso{
  \code{\link{wasserstein}},
  \code{\link{alphaComplexDiag}}, \code{\link{alphaComplexDiag}}, \code{\link{gridDiag}}, \code{\link{ripsDiag}},
  \code{\link{plot.diagram}}
}

\examples{
XX1 <- circleUnif(20)
XX2 <- circleUnif(20, r = 0.2)

DiagLim <- 5
maxdimension <- 1

Diag1 <- ripsDiag(XX1, maxdimension, DiagLim, printProgress = FALSE)
Diag2 <- ripsDiag(XX2, maxdimension, DiagLim, printProgress = FALSE)

bottleneckDist <- bottleneck(Diag1[["diagram"]], Diag2[["diagram"]],
                             dimension = 1)
print(bottleneckDist)
}

\keyword{methods}
\keyword{optimize}
