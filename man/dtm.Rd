\name{dtm}

\alias{dtm}

\title{
  Distance to Measure Function
}

\description{
  The function \code{dtm} computes the "distance to measure function" on a set of points \code{Grid}, using the uniform empirical measure on a set of points \code{X}. Given a probability measure \eqn{P}, The distance to measure function, for each \eqn{y \in R^d}{y in R^d}, is defined by
  \deqn{
    d_{m0}(y) = \left(\frac{1}{m0}\int_0^{m0} ( G_y^{-1}(u))^{r} du\right)^{1/r},
  }{
    d_{m0}(y) = (1/m0 int_0^{m0} (G_y^{-1}(u))^{r} du)^{1/r},
  }
  where \eqn{G_y(t) = P( \Vert X-y \Vert \le t)}{G_y(t) = P(||X-y|| \le t)}, and \eqn{m0 \in (0,1)}{m0 in (0,1)} and \eqn{r \in [1,\infty)}{r in [1,\infty)} are tuning parameters. As \code{m0} increases, DTM function becomes smoother, so \code{m0} can be understood as a smoothing parameter. \code{r} affects less but also changes DTM function as well. The DTM can be seen as a smoothed version of the distance function. See Details and References.

  Given \eqn{X=\{x_1, \dots, x_n\}}{X={x_1, \dots, x_n}}, the empirical version of the distance to measure is
  \deqn{
    \hat d_{m0}(y) = \left(\frac{1}{k} \sum_{x_i \in N_k(y)} \Vert x_i-y \Vert^{r}\right)^{1/r},
  }{
    \hat d_{m0}(y) = (1/k \sum_{x_i in N_k(y)} ||x_i-y||^r)^{1/r},
  }
  where \eqn{k= \lceil m0 * n \rceil} and \eqn{N_k(y)} is the set containing the \eqn{k} nearest neighbors of \eqn{y} among \eqn{x_1, \ldots, x_n}.
}

\usage{
dtm(X, Grid, m0, r = 2, weight = 1)
}

\arguments{
  \item{X}{
    an \eqn{n} by \eqn{d} matrix of coordinates of points used to construct the uniform empirical measure for the distance to measure, where \eqn{n} is the number of points and \eqn{d} is the dimension.
  }

  \item{Grid}{
    an \eqn{m} by \eqn{d} matrix of coordinates of points where the distance to measure is computed, where \eqn{m} is the number of points in \code{Grid} and \eqn{d} is the dimension.
  }

  \item{m0}{
    a numeric variable for the smoothing parameter of the distance to measure. Roughly, \code{m0} is the the percentage of points of \code{X} that are considered when the distance to measure is computed for each point of \code{Grid}. The value of \code{m0} should be in \eqn{(0,1)}.
  }

  \item{r}{
    a numeric variable for the tuning parameter of the distance to measure. The value of \code{r} should be in \eqn{[1,\infty)}, and the default value is \code{2}.
  }

  \item{weight}{
    either a number, or a vector of length \eqn{n}. If it is a number, then same weight is applied to each points of \code{X}. If it is a vector, \code{weight} represents weights of each points of \code{X}. The default value is \code{1}.
  }
}

\details{
  See (Chazal, Cohen-Steiner, and Merigot, 2011, Definition 3.2) and (Chazal, Massart, and Michel, 2015, Equation (2)) for a formal definition of the "distance to measure" function.
}

\value{
  The function \code{dtm} returns a vector of length \eqn{m} (the number of points stored in \code{Grid}) containing the value of the distance to measure function evaluated at each point of \code{Grid}.
}

\references{
  Chazal F, Cohen-Steiner D, Merigot Q (2011). "Geometric inference for probability measures." Foundations of Computational Mathematics 11.6, 733-751.

  Chazal F, Massart P, Michel B (2015). "Rates of convergence for robust geometric inference."

  Chazal F, Fasy BT, Lecci F, Michel B, Rinaldo A, Wasserman L (2014). "Robust Topological Inference: Distance-To-a-Measure and Kernel Distance." Technical Report. 
}

\author{
  Jisu Kim and Fabrizio Lecci
}

\seealso{
  \code{\link{kde}}, \code{\link{kernelDist}}, \code{\link{distFct}}
}

\examples{
## Generate Data from the unit circle
n <- 300
X <- circleUnif(n)

## Construct a grid of points over which we evaluate the function
by <- 0.065
Xseq <- seq(-1.6, 1.6, by = by)
Yseq <- seq(-1.7, 1.7, by = by)
Grid <- expand.grid(Xseq, Yseq)

## distance to measure
m0 <- 0.1
DTM <- dtm(X, Grid, m0)
}

\keyword{nonparametric}
