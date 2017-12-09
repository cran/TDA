kde <-
function(X, Grid, h, kertype = "Gaussian", weight = 1, printProgress = FALSE) {

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.numeric(Grid) && !is.data.frame(Grid)) {
    stop("Grid should be a matrix of coordinates")
  }
  if (NCOL(X) != NCOL(Grid)) {
    stop("dimensions of X and Grid do not match")
  }
  if (!is.numeric(h) || length(h) != 1 || h <= 0) {
    stop("h should be a positive number")
  }
  if (kertype == "gaussian" || kertype == "GAUSSIAN") {
    kertype <- "Gaussian"
  }
  if (kertype != "Gaussian" && kertype != "Epanechnikov") {
    stop("kertype for specifying kernel function should be a string: either 'Gaussian' or 'Epanechnikov'")
  }
  if (!is.numeric(weight) ||
      (length(weight) != 1 && length(weight) != NROW(X))) {
    stop("weight should be either a number or a vector of length equals the number of sample")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be a logical variable")
  }

  return (Kde(X = as.matrix(X), Grid = as.matrix(Grid), h = h,
      kertype = kertype, weight = weight, printProgress = printProgress))
}
