\name{dtm}
\alias{dtm}
\title{
Distance to Measure Function
}
\description{
This function computes the "distance to measure function" on a set of points \code{Grid}, using the uniform empirical measure on a set of points \code{X}. Given a probability measure \eqn{P}, The distance to measure function, for each \eqn{y \in R^d}, is defined by
\deqn{
d_{m_0}(y) =\sqrt{ \frac{1}{m_0}\int_0^{m_0} ( G_y^{-1}(u))^2 du},
}
where
\eqn{G_y(t) = P( \Vert X-y \Vert \leq t)} and \eqn{0<m_0<1} is a smoothing parameter.
See Details and References.

Given \eqn{X=\{x_1, \dots, x_n\}}, the empirical version of the distance to measure is
\deqn{
\hat d_{m_0}(y) = \sqrt{ \frac{1}{k} \sum_{x_i \in N_k(y)} \Vert x_i-y \Vert^2 },
}
where \eqn{k= \lceil m_0 n \rceil} and \eqn{N_k(y)} is the 
set containing the \eqn{k} nearest neighbors 
of \eqn{y} among \eqn{x_1, \ldots, x_n}.
}
\usage{
dtm(X, Grid, m0)
}
\arguments{
  \item{X}{
	an \eqn{n} by \eqn{d} matrix of coordinates of points used to construct the uniform empirical measure for the distance to measure, where \eqn{n} is the number of points and \eqn{d} is the dimension.
}
  \item{Grid}{
an \eqn{m} by \eqn{d} matrix of coordinates, where \eqn{m} is the number of points in \code{Grid}.
}
  \item{m0}{
a numeric variable for the smoothing parameter of the distance to measure. Roughly, \code{m0} is the the percentage of points of \code{X} that are considered when the distance to measure is computed for each point of \code{Grid}.
}
}
\details{
See Definition 3.2 of the reference for a formal definition of the "distance to measure" function.
}
\value{
\code{dtm} returns a vector of length \eqn{m} (the number of points stored in \code{Grid}) containing the value of the distance to measure function evaluated at each point of \code{Grid}.
}
\references{
Frederic Chazal, David Cohen-Steiner, and Quentin Merigot. "Geometric inference for probability measures." Foundations of Computational Mathematics 11.6 (2011): 733-751.

Chazal F, Fasy BT, Lecci F, Michel B, Rinaldo A, Wasserman L (2014). "Robust Topological Inference: Distance-To-a-Measure and Kernel Distance." Technical Report.	
}
\author{
Fabrizio Lecci
}
\seealso{
\code{\link{kde}}, \code{\link{kernelDist}}, \code{\link{distFct}}
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

## distance to measure
m0=0.1
DTM=dtm(X, Grid, m0)
}
\keyword{nonparametric}

