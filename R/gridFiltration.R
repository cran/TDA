gridFiltration <-
function(X = NULL, FUN = NULL, lim = NULL, by = NULL, FUNvalues = NULL,
         maxdimension = max(NCOL(X), length(dim(FUNvalues))) - 1,
         sublevel = TRUE, printProgress = FALSE, ...) {

  if (!xor(is.null(X) || is.null(FUN), is.null(FUNvalues))) {
    stop("either values of X, FUN should be set, or a value of FUNvalues should be set, but not both")
  }
  if (!is.null(X) && !is.null(FUN)) {
    if (!is.numeric(X) && !is.data.frame(X)) {
      stop("X should be a matrix of coordinates")
    }
    if (!is.function(FUN)) {
      stop("FUN should be a function")
    }
    if (is.null(lim) || is.null(by)) {
      stop("When X and FUN are set, lim and by should also be set")
    }
  }
  if (!is.null(lim) && !is.null(by)) {
    if (!is.numeric(lim) || length(lim) %% 2 != 0) {
      stop("lim should be either a numeric matrix or a numeric vector of even length")
    }
    if (!is.numeric(by) || any(by <= 0)) {
      stop("by should be positive")
    }
  }
  if (!is.null(X) && !is.null(FUN) && !is.null(lim) && !is.null(by)) {
    if (2 * NCOL(X) != length(lim)) {
      stop("dimension of X does not match with lim")
    }
    if (length(by) != 1 && length(by) != NCOL(X)) {
      stop("by should be either a number or a vector of length equals dimension of grid")
    }
  }
  if (!is.null(FUNvalues)) {
    if (!is.numeric(FUNvalues) && !is.data.frame(FUNvalues)) {
      stop("FUNvalues should be an array")
    }
  }
  if (!is.numeric(maxdimension) || length(maxdimension) != 1 ||
      maxdimension < 0) {
    stop("maxdimnsion should be a nonnegative integer")
  }
  if (!is.logical(sublevel)) {
    stop("sublevel should be logical")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }

  if (!is.null(X)) {
    X <- as.matrix(X)
  }
  if (!is.null(lim) && !is.null(by)) {
    Grid <- gridBy(lim = lim, by = by)
  }
  if (is.null(FUNvalues)) {
    FUNvalues <- FUN(X, Grid[["grid"]], ...)
    gridDim <- Grid[["dim"]]
  } else {
    if (is.data.frame(FUNvalues)) {
      FUNvalues <- as.matrix(FUNvalues)
    } else {
      FUNvalues <- as.array(FUNvalues)
    }
    gridDim <- dim(FUNvalues)
  }

  if (!is.null(lim) && !is.null(by) &&
      length(Grid[["dim"]]) != length(gridDim)) {
    stop("dimension of FUNvalues does not match with lim and by")
  }

  maxdimension <- min(maxdimension, length(gridDim) - 1)
  if (sublevel == FALSE) {
    FUNvalues <- -FUNvalues
  }

  # compute persistence diagram of function values over a grid
  if (length(gridDim) <= 3) {
    gridOut <- GridFiltration(
	    FUNvalues = FUNvalues, gridDim = as.integer(gridDim),
		maxdimension = as.integer(maxdimension), decomposition = "5tetrahedra",
		printProgress = printProgress)
  } else {
    gridOut <- GridFiltration(
	    FUNvalues = FUNvalues, gridDim = as.integer(gridDim),
        maxdimension = as.integer(maxdimension), decomposition = "barycenter",
		printProgress = printProgress)
  }

  if (sublevel == FALSE) {
    gridOut[[2]] <- -gridOut[[2]]
  }
  
  if (!is.null(lim) && !is.null(by)) {
    out <- list(
	    "cmplx" = gridOut[[1]], "values" = gridOut[[2]],
        "increasing" = sublevel, "coordinates" = Grid[["grid"]])
  } else {
    out <- list(
        "cmplx" = gridOut[[1]], "values" = -gridOut[[2]],
        "increasing" = sublevel)
  }

  return (out)
}