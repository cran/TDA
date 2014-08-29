\name{ripsDiag}
\alias{ripsDiag}
\title{
Rips Persistence Diagram
}
\description{
This function computes the persistence diagram of the Rips filtration built on top of a point cloud.
}
\usage{
ripsDiag(X, maxdimension, maxscale, dist = "euclidean", 
         library="GUDHI", printProgress = FALSE)
}
\arguments{
  \item{X}{
if \code{dist="euclidean"}, \code{X} is an \eqn{n} by \eqn{d} matrix of coordinates, where \eqn{n} is the number of points in the \eqn{d}-dimensional euclidean space.
if \code{dist="arbitrary"}, \code{X} is an \eqn{n} by \eqn{n} matrix of distances of \eqn{n} points.
}
  \item{maxdimension}{
integer: max dimension of the homological features to be computed. (e.g. 0 for connected components, 1 for connected components and loops, 2 for connected components, loops, voids, etc.)
}
  \item{maxscale}{
number: maximum value of the rips filtration.
}
  \item{dist}{
\code{"euclidean"} for Euclidean distance, \code{"arbitrary"} for an arbitrary distance given in input as a distance matrix.
}
  \item{library}{
	if dist=='euclidean', the user can compute the Rips persistence diagram using either the library 'GUDHI', or 'Dionysus'. Default is 'GUDHI'.
}

  \item{printProgress}{
logical: if \code{TRUE}, a progress bar is printed. Default is \code{FALSE}.
}
}
\details{
For Rips Diagrams based on Euclidean distance of the input point cloud, the user can decide to use either the C++ library GUDHI, or the C++ library Dionysus.
For Rips Diagrams based on arbitrary distance, the the C++ library Dionysus is used.
See refereneces.
}
\value{
\code{ripsDiag} returns a list with the following elements:
\item{diagram}{
an object of class \code{diagram}, a \eqn{P} by 3 matrix, where \eqn{P} is the number of points in the resulting persistence diagram. The first column contains the dimension of each feature (0 for components, 1 for loops, 2 for voids, etc.). Second and third columns are Birth and Death of the features.
}
}
\references{
Clement Maria, "GUDHI, Simplicial Complexes and Persistent Homology Packages". https://project.inria.fr/gudhi/software/.

Dmitriy Morozov, "Dionysus, a C++ library for computing persistent homology". http://www.mrzv.org/software/dionysus/

Herbert Edelsbrunner and John Harer (2010), Computational topology: an introduction. American Mathematical Society.

Brittany Fasy, Fabrizio Lecci, Alessandro Rinaldo, Larry Wasserman, Sivaraman Balakrishnan, and Aarti Singh. (2013), "Statistical Inference For Persistent Homology", (arXiv:1303.7117). To appear, Annals of Statistics.
}
\author{
Brittany T. Fasy, Fabrizio Lecci and Clement Maria
}
\seealso{
\code{\link{summary.diagram}}, \code{\link{plot.diagram}},
\code{\link{gridDiag}}
}
\examples{
## EXAMPLE 1: rips diagram for circles (euclidean distance)
XX = circleUnif(30)
maxscale=5
maxdimension=1
## note that the input XX is a point cloud
Diag=ripsDiag(XX,maxdimension,maxscale, printProgress=TRUE)


## EXAMPLE 2: rips diagram with arbitrary distance
## distance matrix for triangle with edges of length: 1,2,4
distX=matrix(c(0,1,2,1,0,4,2,4,0), ncol=3)
maxscale=5
maxdimension=1
## note that the input distXX is a distance matrix
DiagTri=ripsDiag(distX,maxdimension,maxscale, dist="arbitrary", printProgress=TRUE)
#points with lifetime=0 are not shown. e.g. the loop of the triangle.
print(DiagTri$diagram)
}
\keyword{methods}

