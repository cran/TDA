\name{plot.maxPersistence}

\alias{plot.maxPersistence}

\title{
  Summary plot for the maxPersistence function
}

\description{
  This function plots an object of class \code{maxPersistence}, for the selection of the optimal smoothing parameter for persistent homology.
  For each value of the smoothing parameter, the plot shows the number of detected features, their persistence, and a bootstrap confidence band.
}

\usage{
\method{plot}{maxPersistence}(x, features = "dimension", colorBand = "pink", colorBorder = NA, ...)
}

\arguments{
  \item{x}{
    an object of class \code{maxPersistence}, as returned by the functions \code{\link{maxPersistence}}
  }

  \item{features}{
    string: if "all" then all the features are plotted; if "dimension" then only the features of the dimension used to compute the confidence band are plotted.
  }

  \item{colorBand}{
    the color for filling the confidence band. The default is "pink". (NA leaves the band unfilled)
  }

  \item{colorBorder}{
    the color to draw the border of the confidence band. The default is NA and omits the border.
  }

  \item{\dots}{
    additional graphical parameters.
  }
}

\references{
  Frederic Chazal,  Jessi Cisewski, Brittany T. Fasy, Fabrizio Lecci, Bertrand Michel, Alessandro Rinaldo, and Larry Wasserman, (2014), "Robust Topological Inference: distance-to-a-measure and kernel distance"

  Brittany T. Fasy, Fabrizio Lecci, Alessandro Rinaldo, Larry Wasserman, Sivaraman Balakrishnan, and Aarti Singh. (2013), "Statistical Inference For Persistent Homology", (arXiv:1303.7117). To appear, Annals of Statistics.
}

\author{
  Fabrizio Lecci
}

\seealso{
  \code{\link{maxPersistence}}
}

\examples{
## input data: circle with clutter noise
n <- 600
percNoise <- 0.1
XX1 <- circleUnif(n)
noise <- cbind(runif(percNoise * n, -2, 2), runif(percNoise * n, -2, 2))
X <- rbind(XX1, noise)

## limits of the Gird at which the density estimator is evaluated
Xlim <- c(-2, 2)
Ylim <- c(-2, 2)
lim <- cbind(Xlim, Ylim)
by <- 0.2

B <- 80
alpha <- 0.05

## candidates
parametersKDE <- seq(0.1, 0.5, by = 0.2)

maxKDE <- maxPersistence(kde, parametersKDE, X, lim = lim, by = by,
                         bandFUN = "bootstrapBand", B = B, alpha = alpha,
                         parallel = FALSE, printProgress = TRUE)
print(summary(maxKDE))

par(mfrow = c(1, 2))
plot(X, pch = 16, cex = 0.5, main = "Circle")
plot(maxKDE)
}

\keyword{hplot}
