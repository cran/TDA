\name{summary.diagram}
\alias{summary.diagram}
\alias{print.summary.diagram}
\alias{print.diagram}
\title{
\code{print} and \code{summary} for \code{diagram}
}
\description{
\code{print.diagram} prints a persistence diagram, a \eqn{P} by 3 matrix, where \eqn{P} is the number of points in the diagram. The first column contains the dimension of each feature (0 for components, 1 for loops, 2 for voids, etc.). Second and third columns are Birth and Death of the features. 

\code{summary.diagram} produces basic summaries of a persistence diagrams.
}
\usage{
\method{print}{diagram}(x, ...)
\method{summary}{diagram}(object, ...)
}
\arguments{
  \item{x}{
	an object of class \code{diagram}
}
  \item{object}{
	an object of class \code{diagram}
}
  \item{\dots}{
   additional arguments affecting the summary produced.
}
}
\author{
Fabrizio Lecci
}
\seealso{
\code{\link{plot.diagram}}, \code{\link{gridDiag}}, \code{\link{ripsDiag}},
}
\examples{
# Generate data from 2 circles
XX1 = circleUnif(30)
XX2 = circleUnif(30, r=2) +3
XX=rbind(XX1,XX2)

DiagLim=5            # limit of the filtration
maxdimension=1       # computes betti0 and betti1

Diag=ripsDiag(XX,maxdimension,DiagLim, printProgress=TRUE)

print(Diag$diagram)
print(summary(Diag$diagram))
}
