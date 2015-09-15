\name{clusterTree}
\alias{clusterTree}
\alias{print.clusterTree}
\title{
Density clustering: the cluster tree
}
\description{
Given a point cloud, or a matrix of distances, this function computes a density estimator and returns the corresponding cluster tree of superlevel sets (lambda tree and kappa tree; see references).
}
\usage{
clusterTree(X, k, h = NULL, density = "knn", dist = "euclidean", d = NULL,
            Nlambda = 100, printProgress = FALSE)
}
\arguments{
  \item{X}{
    If \code{dist = "euclidean"} then \code{X} is an \eqn{n} by \eqn{d} matrix of coordinates, where \eqn{n} is the number of points stored in \code{X} and \eqn{d} is the dimension of the space.
    If \code{dist = "arbitrary"} then \code{X} is an \eqn{n} by \eqn{n} matrix of distances. Default is \code{"euclidean"}
}
  \item{k}{
    an integer value specifying the parameter of the underlying k-nearest neighbor similarity graph, used to determine connected components. If \code{density = "knn"}, then \code{k} is also used to compute the k-nearest neighbor density estimator.
}
  \item{h}{
    real value: if \code{density = "kde"}, then \code{h} is used to compute the kernel density estimator with bandwidth \code{h}. Default is \code{NULL}.
}
  \item{density}{
    string: if \code{"knn"} then the k-nearest neighbor density estimator is used to compute the cluster tree;
    if \code{"kde"} then the kernel density estimator is used to compute the cluster tree.
    Default is "knn".
}
  \item{dist}{
string: can be \code{"euclidean"}, when \code{X} is a point cloud or \code{"arbitrary"}, when \code{X} is a matrix of distances. Default is 'euclidean'
}
  \item{d}{
    integer: if \code{dist = "arbitrary"}, then \code{d} is the dimension of the underlying space.
}
  \item{Nlambda}{
    integer: size of the grid of values of the density estimator, used to compute the cluster tree. High \code{Nlambda} (i.e. a fine grid) means a more accurate cluster Tree.
}
  \item{printProgress}{
  logical: if TRUE a progress bar is printed. Default is \code{FALSE}.
}
}
\details{
This function is an implementation of Algorithm 1 in the first reference.
}
\value{
This function returns an object of class \code{clusterTree}, a list with the following components
\item{density}{Vector of length \code{n}: the values of the density estimator evaluated at each of the points stored in \code{X}}
\item{DataPoints}{A list whose elements are the points of \code{X} corresponding to each branch, in the same order of \code{id}}
\item{n}{The number of points stored in the input matrix \code{X}}
\item{id}{Vector: the IDs associated to the branches of the cluster tree}
\item{sons}{A list whose elements are the IDs of the sons of each branch, in the same order of \code{id}}
\item{parent}{Vector: the IDs of the parents of each branch, in the same order of \code{id}}
\item{silo}{A list whose elements are the horizontal coordinates of the silo of each branch, in the same order of \code{id}}
\item{Xbase}{Vector: the horiontal coordinates of the branches of the cluster tree, in the same order of \code{id}}
\item{lambdaBottom}{Vector: the vertical bottom coordinates of the branches of the lambda tree, in the same order of \code{id}}
\item{lambdaTop}{Vector: the vertical top coordinates of the branches of the lambda tree, in the same order of \code{id}}
\item{rBottom}{(only if \code{density = "knn"}) Vector: the vertical bottom coordinates of the branches of the r tree, in the same order of \code{id}}
\item{rTop}{(only if \code{density = "knn"}) Vector: the vertical top coordinates of the branches of the r tree, in the same order of \code{id}}
\item{alphaBottom}{Vector: the vertical bottom coordinates of the branches of the alpha tree, in the same order of \code{id}}
\item{alphaTop}{Vector: the vertical top coordinates of the branches of the alpha tree, in the same order of \code{id}}
\item{Kbottom}{Vector: the vertical bottom coordinates of the branches of the kappa tree, in the same order of \code{id}}
\item{Ktop}{Vector: the vertical top coordinates of the branches of the kappa tree, in the same order of \code{id}}
}
\references{
Brian P. Kent, Alessandro Rinaldo, and Timothy Verstynen, (2013), "DeBaCl: A Python Package for Interactive DEnsity-BAsed CLustering."arXiv:1307.8136 

Fabrizio Lecci, Alessandro Rinaldo, and Larry Wasserman, (2014), "Metric Embeddings for Cluster Trees"
}
\author{
Fabrizio Lecci
}

\seealso{
\code{\link{plot.clusterTree}}
}
\examples{
## Generate data: 3 clusters
n <- 1200    #sample size
Neach <- floor(n / 4) 
X1 <- cbind(rnorm(Neach, 1, .8), rnorm(Neach, 5, 0.8))
X2 <- cbind(rnorm(Neach, 3.5, .8), rnorm(Neach, 5, 0.8))
X3 <- cbind(rnorm(Neach, 6, 1), rnorm(Neach, 1, 1))
X <- rbind(X1, X2, X3)

k <- 100     #parameter of knn

## Density clustering using knn and kde
Tree <- clusterTree(X, k, density = "knn")
TreeKDE <- clusterTree(X, k, h = 0.3, density = "kde")

par(mfrow = c(2, 3))
plot(X, pch = 19, cex = 0.6)
# plot lambda trees
plot(Tree, type = "lambda", main = "lambda Tree (knn)")
plot(TreeKDE, type = "lambda", main = "lambda Tree (kde)")
# plot clusters
plot(X, pch = 19, cex = 0.6, main = "cluster labels")
for (i in Tree[["id"]]){
  points(matrix(X[Tree[["DataPoints"]][[i]],],ncol = 2), col = i, pch = 19,
         cex = 0.6)
}
#plot kappa trees
plot(Tree, type = "kappa", main = "kappa Tree (knn)")
plot(TreeKDE, type = "kappa", main = "kappa Tree (kde)")
}
\keyword{nonparametric}
