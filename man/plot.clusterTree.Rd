\name{plot.clusterTree}

\alias{plot.clusterTree}

\title{
  Plots the Cluster Tree
}

\description{
  The function \code{plot.clusterTree} plots the Cluster Tree stored in an object of class \code{clusterTree}.
}

\usage{
\method{plot}{clusterTree}(
    x, type = "lambda", color = NULL, add = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{x}{
    an object of class \code{clusterTree}. (see \code{\link{clusterTree}})
  }

  \item{type}{
    string: if "lambda", then the lambda Tree is plotted.
            if "r", then the r Tree is plotted.
            if "alpha", then the alpha Tree is plotted.
            if "kappa", then the kappa Tree is plotted.
  }

  \item{color}{
    number: the color of the branches of the Cluster Tree. The default value is \code{NULL} and a different color is assigned to each branch.
  }

  \item{add}{
    logical: if \code{TRUE}, the Tree is added to an existing plot.
  }

  \item{\dots}{
    additional graphical parameters.
  }
}

\references{
  Kent BP, Rinaldo A, Verstynen T (2013). "DeBaCl: A Python Package for Interactive DEnsity-BAsed CLustering." arXiv:1307.8136 

  Lecci F, Rinaldo A, Wasserman L (2014). "Metric Embeddings for Cluster Trees"
}

\author{
  Fabrizio Lecci
}

\seealso{
  \code{\link{clusterTree}}, \code{\link{print.clusterTree}}
}

\examples{
## Generate data: 3 clusters
n <- 1200  #sample size
Neach <- floor(n / 4) 
X1 <- cbind(rnorm(Neach, 1, .8), rnorm(Neach, 5, 0.8))
X2 <- cbind(rnorm(Neach, 3.5, .8), rnorm(Neach, 5, 0.8))
X3 <- cbind(rnorm(Neach, 6, 1), rnorm(Neach, 1, 1))
XX <- rbind(X1, X2, X3)

k <- 100   #parameter of knn

## Density clustering using knn and kde
Tree <- clusterTree(XX, k, density = "knn")
TreeKDE <- clusterTree(XX,k, h = 0.3, density = "kde")

par(mfrow = c(2, 3))
plot(XX, pch = 19, cex = 0.6)
# plot lambda trees
plot(Tree, type = "lambda", main = "lambda Tree (knn)")
plot(TreeKDE, type = "lambda", main = "lambda Tree (kde)")
# plot clusters
plot(XX, pch = 19, cex = 0.6, main = "cluster labels")
for (i in Tree[["id"]]){
  points(matrix(XX[Tree[["DataPoints"]][[i]], ], ncol = 2), col = i, pch = 19,
         cex = 0.6)
}
#plot kappa trees
plot(Tree, type = "kappa", main = "kappa Tree (knn)")
plot(TreeKDE, type = "kappa", main = "kappa Tree (kde)")
}

\keyword{hplot}
