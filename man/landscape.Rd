\name{landscape}
\alias{landscape}
\title{
The Persistence Landscape Function
}
\description{
This function computes the landscape function corresponding to a given persistence diagram.
}
\usage{
landscape(Diag, dimension = 1, KK = 1, 
          tseq=seq( min(Diag[,2:3]), max(Diag[,2:3]), length=500) )
}
\arguments{
  \item{Diag}{
	an object of class \code{diagram} or a \eqn{P} by \eqn{3} matrix, storing a persistence diagram with columns: "dimension", "Birth", "Death". 
}
  \item{dimension}{
	the dimension of the topological features under consideration. Default is 1 (loops).
}
  \item{KK}{
	the order of the landscape function. Default is 1. (First Landscape function).
}
  \item{tseq}{
	a vector of values at which the landscape function is evaluated.
}
}
\value{
	Returns a numeric vector of the same length of \code{tseq}, with the values of the landscape function evaluated at each point of \code{tseq}.
}
\references{
Peter Bubenik, (2012), "Statistical topology using persistence landscapes", arXiv1207.6437.


Frederic Chazal, Brittany T. Fasy, Fabrizio Lecci, Alessandro Rinaldo, and Larry Wasserman, (2014), "Stochastic Convergence of Persistence Landscapes and Silhouettes", Proceedings of the 30th Symposium of Computational Geometry (SoCG). (arXiv:1312.0308)
}
\author{
Fabrizio Lecci
}
\seealso{
\code{\link{silhouette}}
}
\examples{
XX1 = circleUnif(30)
XX2 = circleUnif(30, r=2) +3
XX=rbind(XX1,XX2)

DiagLim=5
maxdimension=1
Diag=ripsDiag(XX,maxdimension,DiagLim, printStatus=TRUE)

#persistence landscape
tseq=seq(0,DiagLim, length=1000)
Land= landscape(Diag, dimension=1, KK=1, tseq)

plot(tseq, Land, type="l", xlab="t", ylab="landscape", asp=1)
}
\keyword{methods}
