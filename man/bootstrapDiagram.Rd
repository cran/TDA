\name{bootstrapDiagram}
\alias{bootstrapDiagram}
\title{
Bootstrapped Confidence Set for a Persistence Diagram, using the Bottleneck Distance (or the Wasserstein distance).
}
\description{
\code{bootstrapDiagram} computes a \code{(1-alpha)} confidence set for the Persistence Diagram of a filtration of sublevel sets (or superlevel sets) of a function evaluated over a grid of points. The function returns the (\code{1-alpha}) quantile of \code{B} bottleneck distances (or Wasserstein distances), computed in \code{B} iterations of the bootstrap algorithm. The method is discussed in the 1st reference. 
}
\usage{
bootstrapDiagram(X, FUN, lim, by, sublevel = TRUE, library="Dionysus",
                 B=30, alpha=0.05, distance="bottleneck", dimension=1,
                 p=1, printProgress = FALSE, ...)
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
  \item{sublevel}{
a logical variable indicating if the Persistence Diagram should be computed for sublevel sets (\code{TRUE}) or superlevel sets (\code{FALSE}) of the function. Default is \code{TRUE}.
}
  \item{library}{
  The user can compute the persistence diagram using either the library 'Dionysus', or 'PHAT'. Default is 'Dionysus'.
}
  \item{B}{
  the number of bootstrap iterations.
}
  \item{alpha}{
  \code{bootstrapDiagram} returns a (1-\code{alpha}) quantile.
}
  \item{distance}{
  a string specifying the distance to be used for persistence diagrams: either 'bottleneck' or 'wasserstein'
}
  \item{dimension}{
	if \code{distance=="wasserstein"}, then \code{dimension} an integer specifying the dimension of the features used to compute the bottleneck distance. 0 for connected components, 1 for loops, 2 for voids and so on.
}
  \item{p}{
		if \code{distance=="wasserstein"}, then \code{p} is an integer specifying the power to be used in the computation of the Wasserstein distance. Default is 1.
}
  \item{printProgress}{
	if \code{TRUE} a progress bar is printed. Default is \code{FALSE}.
}
  \item{...}{
	additional parameters for the function \code{FUN}.
}

}
\details{
\code{bootstrapDiagram} uses \code{gridDiag} to compute the persistence diagram of the input function using the entire sample. Then the bootstrap algorithm, for \code{B} times, computes the bottleneck distance between the original persistence diagram and the one computed using a subsample. Finally the (\code{1-alpha}) quantile of these \code{B} values is returned.
}
\value{
Returns the (1-\code{alpha}) quantile of the values computed by the bootstrap algorithm. 
}

\references{
Chazal F, Fasy BT, Lecci F, Michel B, Rinaldo A, Wasserman L (2014). "Robust Topological Inference: Distance-To-a-Measure and Kernel Distance." Technical Report.	

Larry Wasserman (2004), "All of statistics: a concise course in statistical inference", Springer.

Dmitriy Morozov, "Dionysus, a C++ library for computing persistent homology". http://www.mrzv.org/software/dionysus/
}
\author{
Fabrizio Lecci
}
\note{
 This function uses the C++ library Dionysus for the computation of bottleneck and Wasserstein distances. See references.
}

\seealso{
\code{\link{bottleneck}}, \code{\link{bootstrapBand}},
\code{\link{distFct}}, \code{\link{kde}}, \code{\link{kernelDist}}, \code{\link{dtm}},
\code{\link{summary.diagram}}, \code{\link{plot.diagram}},
}
\examples{
## confidence set for the Kernel Density Diagram

# input data
n = 400
XX = circleUnif(n)

## Ranges of the grid
Xlim=c(-1.8,1.8)
Ylim=c(-1.6,1.6)
lim=cbind(Xlim, Ylim)
by=0.05

h = .3  #bandwidth for the function kde

#Kernel Density Diagram of the superlevel sets
Diag=gridDiag(XX, kde, lim=lim, by=by, sublevel=FALSE, printProgress=TRUE, h=h) 

# confidence set
B=5         ## the number of bootstrap iterations should be higher! 
             ## this is just an example
alpha=0.05

cc=bootstrapDiagram(XX, kde, lim=lim, by=by, sublevel=FALSE, B=B, alpha=alpha, 
   dimension=1, printProgress=TRUE, h=h)

plot(Diag$diagram, band=2*cc)
}
\keyword{nonparametric}
\keyword{htest}
