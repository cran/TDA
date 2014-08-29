\name{gridDiag}
\alias{gridDiag}
\title{
Persistence Diagram of a function over a Grid
}
\description{
\code{gridDiag} computes the Persistence Diagram of a filtration of sublevel sets (or superlevel sets) of a function evaluated over a grid of points in arbitrary dimension \code{d}.
}
\usage{
gridDiag(X, FUN, lim, by, maxdimension=length(lim)/2-1, sublevel = TRUE,
         library="Dionysus", printProgress = FALSE, diagLimit=NULL, ...)
}
\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates, used by the function \code{FUN}, where \eqn{n} is the number of points stored in \code{X} and \eqn{d} is the dimension of the space.
  }
  \item{FUN}{
    a function whose inputs are 1) an \eqn{n} by \eqn{d} matrix of coordinates \code{X}, 2) an \eqn{m} by \eqn{d} matrix of coordinates \code{Grid}, 3) an optional smoothing parameter, and returns a numeric vector of length \eqn{m}. For example see \code{\link{distFct}},  \code{\link{kde}}, and \code{\link{dtm}}  which compute the distance function, the kernel density estimator and the distance to measure, over a grid of points using the input \code{X}. Note that \code{Grid} is not an input of \code{gridDiag}, but is automatically computed by the function using \code{lim}, and \code{by}. 
  }  
  \item{lim}{
    a \eqn{2} by \eqn{d} matrix, where each column specifying the range of each dimension of the grid, over which the function \code{FUN} is evaluated.
  }
  \item{by}{
    either a number or a vector of length \eqn{d} specifying space between points of the grid in each dimension. If a number is given, then same space is used in each dimension.
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
The function evaluates the function \code{FUN} over a grid. Then it constructs a filtration by triangulating the grid and considering the simplices determined by the values of the function of dimension up to \code{maxdimension+1}.
}
\value{
\code{gridDiag} returns a list with the following elements:

\item{diagram}{
an object of class \code{diagram}, a \eqn{P} by 3 matrix, where \eqn{P} is the number of points in the resulting persistence diagram. The first column contains the dimension of each feature (0 for components, 1 for loops, 2 for voids, etc). Second and third columns are Birth and Death of features, in case of a filtration constructed using sublevel sets, or Death and Birth of features, in case of a filtration constructed using superlevel sets.
}
}
\references{
Brittany Fasy, Fabrizio Lecci, Alessandro Rinaldo, Larry Wasserman, Sivaraman Balakrishnan, and Aarti Singh. (2013), "Statistical Inference For Persistent Homology", (arXiv:1303.7117). To appear, Annals of Statistics.

Dmitriy Morozov, "Dionysus, a C++ library for computing persistent homology". http://www.mrzv.org/software/dionysus/

Ulrich Bauer, Michael Kerber, Jan Reininghaus, "PHAT, a software library for persistent homology". https://code.google.com/p/phat/
}
\author{
Brittany T. Fasy, Jisu Kim, Fabrizio Lecci
}
\note{
  The user can decide to use either the C++ library Dionysus or the C++ library PHAT. See references.
 Since dimension of simplicial complex from grid points in \eqn{R^d} is up to \eqn{d}, homology of dimension \eqn{>=d} is trivial. Hence setting \code{maxdimension} with values \eqn{>=d} is equivalent to \code{maxdimension=col(X)-1}.
}

\seealso{
\code{\link{summary.diagram}}, \code{\link{plot.diagram}},
\code{\link{distFct}}, \code{\link{kde}}, \code{\link{kernelDist}}, \code{\link{dtm}}, \code{\link{ripsDiag}}
}
\examples{
## Distance Function Diagram and Kernel Density Diagram

# input data
n = 300
XX = circleUnif(n)

## Ranges of the grid
Xlim=c(-1.8,1.8)
Ylim=c(-1.6,1.6)
lim=cbind(Xlim, Ylim)
by=0.05

h = .3  #bandwidth for the function kde

#Distance Function Diagram of the sublevel sets
Diag1=gridDiag(XX,distFct, lim=lim, by=by, sublevel=TRUE, printProgress=TRUE) 

#Kernel Density Diagram of the superlevel sets
Diag2=gridDiag(XX, kde, lim=lim, by=by, sublevel=FALSE, printProgress=TRUE, h=h) 

#plot
par(mfrow=c(1,3))
plot(XX,cex=0.5, pch=19)
title(main="Data")
plot(Diag1$diagram)
title(main="Distance Function Diagram")
plot(Diag2$diagram)
title(main="Density Persistence Diagram")
}
\keyword{methods}
