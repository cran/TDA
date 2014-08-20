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
maxPersistence(FUN, parameters, X, Xlim, Ylim=NA, Zlim=NA, by, sublevel = TRUE, 
             B = 30, alpha = 0.05, parallel = FALSE, printProgress = FALSE)
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
a logical variable indicating if the persistent homology should be computed for sublevel sets of \code{FUN} (TRUE) or superlevel sets (FALSE). Default is TRUE.

}
  \item{B}{
  the number of bootstrap iterations.
}
  \item{alpha}{
  for each value store in \code{parameters}, \code{maxPersistence} computes a (1-\code{alpha}) confidence band.
}
  \item{parallel}{
    logical: if \code{TRUE} the bootstrap iterations are parallelized, using the library \code{parallel}.
}
  \item{printProgress}{
	if TRUE a progress bar is printed. Default is \code{FALSE}.
}
}
\details{
	\code{maxPersistence} calls the \code{\link{gridDiag}} function, which computes the persistence diagram of sublevel (or superlevel) sets of a function, evaluated over a grid of points in dimension 1,2, or 3.
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
by=0.2

B=80
alpha=0.05

## candidates
parametersKDE=seq(0.1,0.5, by=0.2)

maxKDE=maxPersistence(kde, parametersKDE, X, Xlim, Ylim, Zlim=NA, by=by, B=B,
                  alpha=alpha, parallel=FALSE, printProgress = TRUE)
print(summary(maxKDE))

par(mfrow=c(1,2))
plot(X, pch=16, cex=0.5, main="Circle")
plot(maxKDE)
}
\keyword{methods}