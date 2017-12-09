\name{multipBootstrap}

\alias{multipBootstrap}

\title{
  Multiplier Bootstrap for Persistence Landscapes and Silhouettes
}

\description{
  The function \code{multipBootstrap} computes a confidence band for the average landscape (or the average silhouette) using the multiplier bootstrap.
}

\usage{
multipBootstrap(
    Y, B = 30, alpha = 0.05, parallel = FALSE,
    printProgress = FALSE)
}

\arguments{
  \item{Y}{
    an \eqn{N} by \eqn{m} matrix of values of \eqn{N} persistence landscapes (or silhouettes) evaluated over a 1 dimensional grid of length \eqn{m}.
  }

  \item{B}{
    the number of bootstrap iterations.
  }

  \item{alpha}{
    \code{multipBootstrap} returns a 1-\code{alpha} confidence band for the mean landscape (or silhouette).
  }

  \item{parallel}{
    logical: if \code{TRUE} the bootstrap iterations are parallelized, using the library \code{parallel}.
  }

  \item{printProgress}{
    logical: if \code{TRUE} a progress bar is printed. The default value is \code{FALSE}.
  }
}
\details{
  See Algorithm 1 in the reference.
}

\value{
  The function \code{multipBootstrap} returns a list with the following elements:

  \item{width}{
    number: half of the width of the unfiorm confidence band; that is, the distance of the upper and lower limits of the band from the empirical average landscape (or silhouette).
  }

  \item{mean}{
    a numeric vector of length \eqn{m}, storing the values of the empirical average landscape (or silhouette) over a 1 dimensional grid of length \eqn{m}.
  }

  \item{band}{
    an \eqn{m} by 2 matrix that stores the values of the lower limit of the confidence band (first column) and upper limit of the confidence band (second column), evaluated over a 1 dimensional grid of length \eqn{m}.
  }
}

\references{
  Chazal F, Fasy BT, Lecci F, Rinaldo A, Wasserman L (2014). "Stochastic Convergence of Persistence Landscapes and Silhouettes." Proceedings of the 30th Symposium of Computational Geometry (SoCG). (arXiv:1312.0308)
}

\author{
Fabrizio Lecci
}

\seealso{
  \code{\link{landscape}}, \code{\link{silhouette}}
}

\examples{
nn <- 3000  #large sample size
mm <- 50    #small subsample size
NN <- 5     #we will compute NN diagrams using subsamples of size mm

XX <- circleUnif(nn)  ## large sample from the unit circle

DiagLim <- 2
maxdimension <- 1
tseq <- seq(0, DiagLim, length = 1000)

Diags <- list()  #here we will store the NN rips diagrams
                 #constructed using different subsamples of mm points
#here we'll store the landscapes
Lands <- matrix(0, nrow = NN, ncol = length(tseq))

for (i in seq_len(NN)){
  subXX <- XX[sample(seq_len(nn), mm), ]
  Diags[[i]] <- ripsDiag(subXX, maxdimension, DiagLim)
  Lands[i, ] <- landscape(Diags[[i]][["diagram"]], dimension = 1, KK = 1, tseq)
}

## now we use the NN landscapes to construct a confidence band
B <- 50
alpha <- 0.05
boot <- multipBootstrap(Lands, B, alpha)

LOWband <- boot[["band"]][, 1]
UPband <- boot[["band"]][, 2]
MeanLand <- boot[["mean"]]

plot(tseq, MeanLand, type = "l", lwd = 2, xlab = "", ylab = "",
     main = "Mean Landscape with band", ylim = c(0, 1.2))
polygon(c(tseq, rev(tseq)), c(LOWband, rev(UPband)), col = "pink")
lines(tseq, MeanLand, lwd = 1, col = 2)
}

\keyword{nonparametric}
\keyword{htest}