\name{gridDiag}
\alias{gridDiag}
\title{
Persistence Diagram of a function over a Grid
}
\description{
\code{gridDiag} computes the Persistence Diagram of a filtration of sublevel sets (or superlevel sets) of a function evaluated over a grid of points in dimension \eqn{d=}1, 2 or 3.
}
\usage{
gridDiag(X, FUN, Xlim, Ylim = NA, Zlim = NA, by=(Xlim[2]-Xlim[1])/20, 
         sublevel = TRUE, printStatus = FALSE, diagLimit=NULL, ...)
}
\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates, used by the function \code{FUN}, where \eqn{n} is the number of points stored in \code{X} and \eqn{d} is the dimension (1, 2 or 3).
  }
  \item{FUN}{
    a function whose inputs are 1) an \eqn{n} by \eqn{d} matrix of coordinates \code{X}, 2) an \eqn{m} by \eqn{d} matrix of coordinates \code{Grid}, 3) an optional smoothing parameter, and returns a numeric vector of length \eqn{m}. For example see \code{\link{distFct}},  \code{\link{kde}}, and \code{\link{dtm}}  which compute the distance function, the kernel density estimator and the distance to measure, over a grid of points using the input \code{X}. Note that \code{Grid} is not an input of \code{gridDiag}, but is automatically computed by the function using \code{Xlim}, \code{Ylim}, \code{Zlim}, and \code{by}. 
  }  
  \item{Xlim}{
    a numeric vector of length 2, specifying the range of the first dimension of the grid, over which the function \code{FUN} is evaluated.
  }

  \item{Ylim}{
    a numeric vector of length 2, specifying the range of the second dimension of the grid, over which the function \code{FUN} is evaluated. \code{NA} for a 1 dimensional grid.
  }

  \item{Zlim}{
    a numeric vector of length 2, specifying the range of the third dimension of the grid, over which the function \code{FUN} is evaluated. \code{NA} for a 1 dimensional or 2 dimensional grid.
  }

  \item{by}{
    number: space between points of the grid in each dimension.    
  }

  \item{sublevel}{
a logical variable indicating if the Persistence Diagram should be computed for sublevel sets (\code{TRUE}) or superlevel sets (\code{FALSE}) of the function. Default is \code{TRUE}.
}
  \item{printStatus}{
	if \code{TRUE} a progress bar is printed. Default is \code{FALSE}.
}

  \item{diagLimit}{
	a number that replaces \code{Inf} (if \code{sublevel} is \code{TRUE}) or \code{-Inf} (if \code{sublevel} is \code{FALSE}) in the Death value of the most persistent connected component. Deafult is \code{NULL} and the min/max of the function is used.
}

  \item{...}{
	additional parameters for the function \code{FUN}.
}

}
\details{
The function evaluates the function \code{FUN} over a grid. Then it constructs a filtration by triangulating the grid and considering the simplices determined by the values of the function.
}
\value{
\code{gridDiag} returns an object of class \code{diagram}, a \eqn{P} by 3 matrix, where \eqn{P} is the number of points in the resulting persistence diagram. The first column contains the dimension of each feature (0 for components, 1 for loops, 2 for voids). Second and third columns are Birth and Death of features, in case of a filtration constructed using sublevel sets, or Death and Birth of features, in case of a filtration constructed using superlevel sets.
}
\references{
Brittany Fasy, Fabrizio Lecci, Alessandro Rinaldo, Larry Wasserman, Sivaraman Balakrishnan, and Aarti Singh. (2013), "Statistical Inference For Persistent Homology", (arXiv:1303.7117). To appear, Annals of Statistics.

http://www.mrzv.org/software/dionysus/
}
\author{
Brittany T. Fasy and Fabrizio Lecci
}
\note{
 This function uses the C++ library Dionysus. See references.
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
by=0.05

h = .3  #bandwidth for the function kde

#Distance Function Diagram of the sublevel sets
Diag1=gridDiag(XX,distFct, Xlim, Ylim, by=by, sublevel=TRUE, printStatus=TRUE) 

#Kernel Density Diagram of the superlevel sets
Diag2=gridDiag(XX, kde, Xlim, Ylim, by=by, sublevel=FALSE, printStatus=TRUE, h=h) 

#plot
par(mfrow=c(1,3))
plot(XX,cex=0.5, pch=19)
title(main="Data")
plot(Diag1)
title(main="Distance Function Diagram")
plot(Diag2)
title(main="Density Persistence Diagram")
}
\keyword{methods}
