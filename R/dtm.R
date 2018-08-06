dtm <-
function(X, Grid, m0, r = 2, weight = 1) {

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.numeric(Grid) && !is.data.frame(Grid)) {
    stop("Grid should be a matrix of coordinates")
  }
  if (NCOL(X) != NCOL(Grid)) {
    stop("dimensions of X and Grid do not match")
  }
  if (!is.numeric(m0) || length(m0) != 1 || m0 < 0 || m0 > 1) {
    stop("m0 should be a number between 0 and 1")
  }
  if (!is.numeric(r) || length(r) != 1 || r < 1) {
    stop("r should be a number greater than or equal to 1")
  }
  if (!is.numeric(weight) || 
      (length(weight) != 1 && length(weight) != NROW(X))) {
    stop("weight should be either a number or a vector of length equals the number of sample")
  }

  # without weight
  if (length(weight) == 1) {
    X <- as.matrix(X)
    weightBound <- m0 * NROW(X)
    knnDistance <- FNN::knnx.dist(
		data = X, query = as.matrix(Grid), k = ceiling(weightBound),
		algorithm = c("kd_tree"))
    return (Dtm(knnDistance = knnDistance, weightBound = weightBound, r = r))

  # with weight
  } else {
    X0 <- as.matrix(X[weight != 0, , drop = FALSE]) 
    weight0 <- weight[weight != 0]
    weight0sort <- sort(weight0)
    weightBound <- m0 * sum(weight0)
    weightSumTemp <- 0
    for (k0 in seq(along = weight0)) {
      weightSumTemp <- weightSumTemp + weight0sort[k0]
      if (weightSumTemp >= weightBound) {
        break
      }
    }
    knnDistanceIndex <- FNN::get.knnx(
	    data = X0, query = as.matrix(Grid), k = k0, algorithm = c("kd_tree"))
    return (DtmWeight(
	    knnDistance = knnDistanceIndex[["nn.dist"]], weightBound = weightBound,
		r = r, knnIndex = knnDistanceIndex[["nn.index"]], weight = weight0))
  }
}