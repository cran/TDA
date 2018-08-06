\name{plot.diagram}

\alias{plot.diagram}

\title{
  Plot the Persistence Diagram
}

\description{
  The function \code{plot.diagram} plots the Persistence Diagram stored in an object of class \code{diagram}. Optionally, it can also represent the diagram as a persistence barcode.
}

\usage{
\method{plot}{diagram}(
    x, diagLim = NULL, dimension = NULL, col = NULL,
    rotated = FALSE, barcode = FALSE, band = NULL, lab.line = 2.2,
    colorBand = "pink", colorBorder = NA, add = FALSE, ...)
}

\arguments{
  \item{x}{
    an object of class \code{diagram} (as returned by the functions \code{\link{alphaComplexDiag}}, \code{\link{alphaComplexDiag}}, \code{\link{gridDiag}}, or \code{\link{ripsDiag}}) or an \eqn{n} by 3 matrix, where \eqn{n} is the number of features to be plotted.
  }

  \item{diagLim}{
    numeric vector of length 2, specifying the limits of the plot. If \code{NULL} then it is automatically computed using the lifetimes of the features.
  }

  \item{dimension}{
    number specifying the dimension of the features to be plotted. If \code{NULL} all the features are plotted.
  }

  \item{col}{
    an optional vector of length \eqn{P} that stores the colors of the topological features to be plotted, where \eqn{P} is the number of topological features stored in \code{x}.
  }

  \item{rotated}{
    logical: if \code{FALSE} the plotted diagram has axes (birth, death), if \code{TRUE} the plotted diagram has axes ((birth+death)/2,(death-birth)/2). The default value is \code{FALSE}. 
  }

  \item{barcode}{
    logical: if \code{TRUE} the persistence barcode is plotted, in place of the diagram.
  }

  \item{band}{
    numeric: if \code{band}!=\code{NULL}, a pink band of size \code{band} is added around the diagonal. If also \code{barcode} is \code{TRUE}, then bars shorter than \code{band} are dotted. The default value is \code{NULL}.
  }

  \item{lab.line}{
    number of lines from the plot edge, where the labels will be placed. The default value is \code{2.2}.
  }

  \item{colorBand}{
    the color for filling the confidence band. The default value is \code{"pink"}. (NA leaves the band unfilled)
  }

  \item{colorBorder}{
    the color to draw the border of the confidence band. The default value is \code{NA} and omits the border.
  }

  \item{add}{
    logical: if \code{TRUE}, the points of \code{x} are added to an existing plot.
  }

  \item{\dots}{
    additional graphical parameters.
  }
}

\references{
  Brittany T. Fasy, Fabrizio Lecci, Alessandro Rinaldo, Larry Wasserman, Sivaraman Balakrishnan, and Aarti Singh. (2013), "Statistical Inference For Persistent Homology", (arXiv:1303.7117). To appear, Annals of Statistics.

  Frederic Chazal, Brittany T. Fasy, Fabrizio Lecci, Alessandro Rinaldo, and Larry Wasserman, (2014), "Stochastic Convergence of Persistence Landscapes and Silhouettes", Proceedings of the 30th Symposium of Computational Geometry (SoCG). (arXiv:1312.0308)
}

\author{
  Fabrizio Lecci
}

\seealso{
  \code{\link{alphaComplexDiag}}, \code{\link{alphaComplexDiag}}, \code{\link{gridDiag}}, \code{\link{ripsDiag}}
}

\examples{
XX1 <- circleUnif(30)
XX2 <- circleUnif(30, r = 2) + 3
XX <- rbind(XX1, XX2)

DiagLim <- 5
maxdimension <- 1

## rips diagram
Diag <- ripsDiag(XX, maxdimension, DiagLim, printProgress = TRUE)

#plot
par(mfrow = c(1, 3))
plot(Diag[["diagram"]])
plot(Diag[["diagram"]], rotated = TRUE)
plot(Diag[["diagram"]], barcode = TRUE)
}

\keyword{hplot}
