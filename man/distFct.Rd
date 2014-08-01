\name{distFct}
\alias{distFct}
\title{
Distance function
}
\description{
This function computes the distance between each point of a set \code{Grid} and the corresponding closest point of another set \code{X}.
}
\usage{
distFct(X, Grid)
}
\arguments{
  \item{X}{
	a numeric \eqn{m} by \eqn{d} matrix of coordinates in the space, where \eqn{m} is the number of points in \code{X} and \eqn{d} is the dimension of the space.
}
  \item{Grid}{
	a numeric \eqn{n} by \eqn{d} matrix of coordinates in the space, where \eqn{n} is the number of points in \code{Grid} and \eqn{d} is the dimension of the space.
}
}
\details{
	Given a set of points \code{X}, the distance function computed at \eqn{g} is defined
	as
	\deqn{
		d(g)= \inf_{x \in X} \| x-g \|_2
	}		
}	
\value{
Returns a numeric vector of length \eqn{n}, where \eqn{n} is the number of points stored in \code{Grid}.
}
\author{
Fabrizio Lecci
}
\seealso{
\code{\link{kde}},\code{\link{kernelDist}}, \code{\link{dtm}}
}
\examples{
## Generate Data from the unit circle
n = 300
X = circleUnif(n)

## Construct a grid of points over which we evaluate the function
by=0.065
Xseq=seq(-1.6, 1.6, by=by)
Yseq=seq(-1.7, 1.7, by=by)
Grid=expand.grid(Xseq,Yseq)

## distance fct
distance= distFct(X, Grid)
}
\keyword{nonparametric}

