\name{funFiltration}

\alias{funFiltration}

\title{
  Filtration from function values
}

\description{
  The function \code{funFiltration} computes the filtration from the complex and the function values.
}

\usage{
funFiltration(FUNvalues, cmplx, sublevel = TRUE)
}

\arguments{
  \item{FUNvalues}{
    The function values on the vertices of the complex.
  }

  \item{cmplx}{
    the complex.
  }

  \item{sublevel}{
    a logical variable indicating if the Persistence Diagram should be computed for sublevel sets (\code{TRUE}) or superlevel sets (\code{FALSE}) of the function. The default value is \code{TRUE}.
  }
}

\details{
  See references.
}

\value{
  The function \code{funFiltration} returns a list with the following elements:

  \item{cmplx}{
    a list representing the complex. Its i-th element represents the vertices of i-th simplex.
  }

  \item{values}{
    a vector representing the filtration values. Its i-th element represents the filtration value of i-th simplex.
  }

  \item{increasing}{
    a logical variable indicating if the filtration values are in increasing order (\code{TRUE}) or in decreasing order (\code{FALSE}).
  }
}

\references{
  Edelsbrunner H, Harer J (2010). "Computational topology: an introduction." American Mathematical Society.
}

\author{
  Jisu Kim
}

\seealso{
  \code{\link{filtrationDiag}}
}

\examples{
n <- 5
X <- cbind(cos(2*pi*seq_len(n)/n), sin(2*pi*seq_len(n)/n))
maxdimension <- 1
maxscale <- 1.5
dist <- "euclidean"
library <- "Dionysus"

FltRips <- ripsFiltration(X = X, maxdimension = maxdimension,
               maxscale = maxscale, dist = "euclidean", library = "Dionysus",
               printProgress = TRUE)

FUNvalues <- X[, 1] + X[, 2]

FltFun <- funFiltration(FUNvalues = FUNvalues, cmplx = FltRips[["cmplx"]])
}
\keyword{methods}
