\name{maxPersistence}
\alias{maxPersistence}
\alias{print.maxPersistence}
\alias{summary.maxPersistence}
\alias{print.summary.maxPersistence}
\title{
Maximal Persistence Method
}
\description{
Given a point cloud and a function built on top of the data, we are interested in studying the evolution of the sublevel sets (or superlevel sets) of the function, using persistent homology. The Maximal Persistence Method selects the optimal smoothing parameter of the function, by maximizing the number of significant topological features, or by maximizing the total significant persistence of the features. For each value of the smoothing parameter, this function computes a persistence diagram using \code{gridDiag} and returns the values of the two criteria, the dimension of detected features, their persistence, and a bootstrapped confidence band. The features that fall outside of the band are statistically significant. See References. 
}
\usage{
maxPersistence(FUN, parameters, X, lim, by, maxdimension=length(lim)/2-1, 
             sublevel = TRUE, library="Dionysus", B=30, alpha=0.05, 
             bandFUN="bootstrapBand", distance="bottleneck", dimension=1, p=1,
             parallel=FALSE, printProgress=FALSE)
}
\arguments{
  \item{FUN}{
    the name of a function whose inputs are: 1) \code{X}, a \eqn{n} by \eqn{d} matrix of coordinates of the input point cloud, where \eqn{d} is the dimension of the space; 2) a matrix of coordinates of points forming a grid at which the function can be evaluated (note that this grid is not passed as an input, but is automatically computed by \code{maxPersistence}); 3) a real valued smoothing parameter. For example, see \code{\link{kde}}, \code{\link{dtm}}, \code{\link{kernelDist}}.
}
  \item{parameters}{
    a numerical vector, storing a sequence of values for the smoothing paramter of \code{FUN} among which \code{maxPersistence} will select the optimal ones.
}
  \item{X}{
    a \eqn{n} by \eqn{d} matrix of coordinates of the input point cloud, where \eqn{d} is the dimension of the space.
}
  \item{lim}{
    a \eqn{2} by \eqn{d} matrix, where each column specifying the range of each dimension of the grid, over which the function \code{FUN} is evaluated.
  }
  \item{by}{
    either a number or a vector of length \eqn{d} specifying space between points of the grid in each dimension. If a number is given, then same space is used in each dimension.
  }
  \item{maxdimension}{
    a number that indicates the maximum dimension to compute persistent homology to. Default is \eqn{d-1}, which is (dimension of embedding space - 1).
  }
  \item{sublevel}{
a logical variable indicating if the persistent homology should be computed for sublevel sets of \code{FUN} (TRUE) or superlevel sets (FALSE). Default is TRUE.
}
  \item{library}{
  User can compute the persistence diagram using either the library 'Dionysus', or 'phat'. Default is 'Dionysus'.
}
  \item{bandFUN}{
  the function to be used in the computation of the confidence band. Either 'bootstrapDiagram' or 'bootstrapBand'.
}

  \item{B}{
  the number of bootstrap iterations.
}
  \item{alpha}{
  for each value store in \code{parameters}, \code{maxPersistence} computes a (1-\code{alpha}) confidence band.
}
  \item{distance}{
  optional (if bandFUN==bootstrapDiagram): a string specifying the distance to be used for persistence diagrams: either 'bottleneck' or 'wasserstein'
}
  \item{dimension}{
	optional (if bandFUN==bootstrapDiagram): an integer specifying the dimension of the features used to compute the bottleneck distance. 0 for connected components, 1 for loops, 2 for voids. Deafault is 1.
}
  \item{p}{
  	optional (if bandFUN==bootstrapDiagram AND distance=='wasserstein'): integer specifying the power to be used in the computation of the Wasserstein distance. Default is 1.
}
  \item{parallel}{
    logical: if \code{TRUE} the bootstrap iterations are parallelized, using the library \code{parallel}. (only if bandFUN=="bootstrapBand")
}
  \item{printProgress}{
	if TRUE a progress bar is printed. Default is \code{FALSE}.
}
}
\details{
	\code{maxPersistence} calls the \code{\link{gridDiag}} function, which computes the persistence diagram of sublevel (or superlevel) sets of a function, evaluated over a grid of points.
}
\value{
  The function returns an object of the class "maxPersistence", a list with the following components
  \item{parameters}{the same vector \code{parameters} given in input
}
  \item{sigNumber}{a numeric vector storing the number of significant features in the persistence diagrams computed using each value in \code{parameters}}
  \item{sigPersistence}{a numeric vector storing the sum of significant persistence of the features in the persistence diagrams, computed using each value in \code{parameters}}
  \item{bands}{a numeric vector storing the bootstrap band's width, 
  	for each value in \code{parameters}}
  \item{Persistence}{a list of the same lenght of \code{parameters}. Each element of the list is a \eqn{P_i} by 2 matrix, where \eqn{P_i} is the number of features
  	found using the parameter \eqn{i}: the first column stores the dimension of each feature and the second column the persistence abs(death-birth|).}

}
\references{
Frederic Chazal,  Jessi Cisewski, Brittany T. Fasy, Fabrizio Lecci, Bertrand Michel, Alessandro Rinaldo, and Larry Wasserman, (2014), "Robust Topological Inference: distance-to-a-measure and kernel distance"

Brittany T. Fasy, Fabrizio Lecci, Alessandro Rinaldo, Larry Wasserman, Sivaraman Balakrishnan, and Aarti Singh. (2013), "Statistical Inference For Persistent Homology", (arXiv:1303.7117). To appear, Annals of Statistics.
}
\author{
Fabrizio Lecci
}
\seealso{
\code{\link{gridDiag}}, \code{\link{kde}}, \code{\link{kernelDist}}, \code{\link{dtm}}, \code{\link{bootstrapBand}}
}
\examples{
## input data: circle with clutter noise
n=600
percNoise=0.1
XX1 = circleUnif(n)
noise=cbind(runif(percNoise*n, -2,2),runif(percNoise*n, -2,2))
X=rbind(XX1,noise)

## limits of the Gird at which the density estimator is evaluated
Xlim=c(-2,2)
Ylim=c(-2,2)
lim=cbind(Xlim,Ylim)
by=0.2

B=80
alpha=0.05

## candidates
parametersKDE=seq(0.1,0.5, by=0.2)

maxKDE=maxPersistence(kde, parametersKDE, X, lim=lim, by=by, 
                  bandFUN="bootstrapBand",B=B, alpha=alpha,
                  parallel=FALSE, printProgress = TRUE)
print(summary(maxKDE))

par(mfrow=c(1,2))
plot(X, pch=16, cex=0.5, main="Circle")
plot(maxKDE)
}
\keyword{methods}