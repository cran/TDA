\name{landscape}

\alias{landscape}

\title{
  The Persistence Landscape Function
}

\description{
  The function \code{landscape} computes the landscape function corresponding to a given persistence diagram.
}

\usage{
landscape(
    Diag, dimension = 1, KK = 1,
    tseq = seq(min(Diag[,2:3]), max(Diag[,2:3]), length=500))
}

\arguments{
  \item{Diag}{
    an object of class \code{diagram} or a \eqn{P} by \eqn{3} matrix, storing a persistence diagram with colnames: "dimension", "Birth", "Death". 
  }

  \item{dimension}{
    the dimension of the topological features under consideration. The default value is \code{1} (loops).
  }

  \item{KK}{
    a vector: the order of the landscape function. The default value is \code{1}. (First Landscape function).
  }

  \item{tseq}{
    a vector of values at which the landscape function is evaluated.
  }
}

\value{
  The function \code{landscape} returns a numeric matrix with the number of row as the length of \code{tseq} and the number of column as the length of \code{KK}. The value at ith row and jth column represents the value of the \code{KK[j]}-th landscape function evaluated at \code{tseq[i]}.
}

\references{
  Bubenik P (2012). "Statistical topology using persistence landscapes." arXiv:1207.6437.

  Chazal F, Fasy BT, Lecci F, Rinaldo A, Wasserman L (2014). "Stochastic Convergence of Persistence Landscapes and Silhouettes." Proceedings of the 30th Symposium of Computational Geometry (SoCG). (arXiv:1312.0308)
}

\author{
  Fabrizio Lecci
}

\seealso{
  \code{\link{silhouette}}
}

\examples{
Diag <- matrix(c(0, 0, 10, 1, 0, 3, 1, 3, 8), ncol = 3, byrow = TRUE)
DiagLim <- 10
colnames(Diag) <- c("dimension", "Birth", "Death")

#persistence landscape
tseq <- seq(0,DiagLim, length = 1000)
Land <- landscape(Diag, dimension = 1, KK = 1, tseq)

par(mfrow = c(1,2))
plot.diagram(Diag)
plot(tseq, Land, type = "l", xlab = "t", ylab = "landscape", asp = 1)
}

\keyword{methods}
