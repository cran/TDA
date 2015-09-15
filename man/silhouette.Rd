\name{silhouette}

\alias{silhouette}

\title{
  The Persistence Silhouette Function
}

\description{
  This function computes the silhouette function corresponding to a given persistence diagram.
}

\usage{
silhouette(Diag, p = 1, dimension = 1, 
           tseq = seq(min(Diag[, 2:3]), max(Diag[, 2:3]), length = 500))
}

\arguments{
  \item{Diag}{
    an object of class \code{diagram} or a \eqn{P} by \eqn{3} matrix, storing a persistence diagram with colnames: "dimension", "Birth", "Death". 
  }

  \item{p}{
    a vector: the power of the weights of the silhouette function. See the definition of silhouette function, Section 5 in the reference.
  }

  \item{dimension}{
    the dimension of the topological features under consideration. Default is 1 (loops).
  }

  \item{tseq}{
    a vector of values at which the silhouette function is evaluated.
  }
}

\value{
  Returns a numeric matrix of with the number of row as the length of \code{tseq} and the number of column as the length of \code{p}. The value at ith row and jth column represents the value of the \code{p[j]}-th power silhouette function evaluated at \code{tseq[i]}.
}

\references{
  Frederic Chazal, Brittany T. Fasy, Fabrizio Lecci, Alessandro Rinaldo, and Larry Wasserman, (2014), "Stochastic Convergence of Persistence Landscapes and Silhouettes", Proceedings of the 30th Symposium of Computational Geometry (SoCG). (arXiv:1312.0308)
}

\author{
  Fabrizio Lecci
}

\seealso{
  \code{\link{landscape}}
}

\examples{
Diag <- matrix(c(0, 0, 10, 1, 0, 3, 1, 3, 8), ncol = 3, byrow = TRUE)
DiagLim <- 10
colnames(Diag) <- c("dimension", "Birth", "Death")

#persistence silhouette
tseq <- seq(0, DiagLim, length = 1000)
Sil <- silhouette(Diag, p = 1,  dimension = 1, tseq)

par(mfrow = c(1, 2))
plot.diagram(Diag)
plot(tseq, Sil, type = "l", xlab = "t", ylab = "silhouette", asp = 1)
}
\keyword{methods}
