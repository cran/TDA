\name{hausdInterval}
\alias{hausdInterval}
\title{
  Subsampling Confidence Interval for the Hausdorff Distance between a Manifold and a Sample
}
\description{
\code{hausdInterval} computes a confidence interval for the Hausdorff distance between a point cloud \code{X} and the underlying manifold from which \code{X} was sampled. See Details. The validity of the method is proved in the 1st Reference.
}
\usage{
hausdInterval(X, m, B = 30, alpha = 0.05, parallel = FALSE, 
              printProgress=FALSE)
}
\arguments{
  \item{X}{
	an \eqn{n} by \eqn{d} matrix of coordinates of sampled points.
}
  \item{m}{
  the size of the subsamples.
}
  \item{B}{
  the number of subsampling iterations.
}
  \item{alpha}{
  \code{hausdInterval} returns a (1-\code{alpha}) confidence interval.
}
  \item{parallel}{
    logical: if \code{TRUE} the iterations are parallelized, using the library \code{parallel}.
}
  \item{printProgress}{
	if TRUE a progress bar is printed. Default is \code{FALSE}.
}
}
\details{
For \code{B} times, the subsampling algorithm subsamples \code{m} points of \code{X} (without replacement) and computes the Hausdorff distance between the original sample \code{X} and the subsample. The result is a sequence of \code{B} values. Let \eqn{q} be the (1-\code{alpha}) quantile of these values and let \eqn{c=2*q}. The interval \eqn{[0, c]} is a valid  (1-\code{alpha}) confidence interval for the Hausdorff distance between \code{X} and the underlying manifold, as proven in Theorem 3 of the first reference.
}
\value{ 
	Returns a number \eqn{c}. The confidence interval is \eqn{[0, c]}. 
}
\references{
Brittany T. Fasy, Fabrizio Lecci, Alessandro Rinaldo, Larry Wasserman, Sivaraman Balakrishnan, and Aarti Singh. (2013), "Statistical Inference For Persistent Homology: Confidence Sets for Persistence Diagrams", (arXiv:1303.7117). To appear, Annals of Statistics.
}
\author{
Fabrizio Lecci
}
\seealso{
\code{\link{bootstrapBand}}
}
\examples{
X= circleUnif(1000)
interval= hausdInterval(X, m=800)
print(interval)
}
\keyword{nonparametric}
\keyword{htest}
