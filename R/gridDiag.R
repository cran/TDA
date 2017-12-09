gridDiag <-
function(X = NULL, FUN = NULL, lim = NULL, by = NULL, FUNvalues = NULL,
         maxdimension = max(NCOL(X), length(dim(FUNvalues))) - 1,
         sublevel = TRUE, library = "GUDHI", location = FALSE,
         printProgress = FALSE, diagLimit = NULL, ...) {

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
  if (library == "gudhi" || library == "Gudhi") {
    library <- "GUDHI"
  }
  if (library == "dionysus" || library == "DIONYSUS") {
    library <- "Dionysus"
  }
  if (library == "phat" || library == "Phat") {
    library <- "PHAT"
  }
  if (library != "GUDHI" && library != "Dionysus" && library != "PHAT") {
    stop("library for computing persistence diagram should be a string: either 'GUDHI', 'Dionysus', or 'PHAT'")
  }
  if (!is.logical(location)) {
    stop("location should be logical")
  }
  if(location == TRUE && (is.null(lim) || is.null(by))) {
    stop("If location is TRUE, both lim and by should be set")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }
  if ((!is.numeric(diagLimit) || length(diagLimit) != 1) &&
      !is.null(diagLimit)) {
    stop("diagLimit should be a positive number")
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
    gridOut <- GridDiag(
	    FUNvalues = FUNvalues, gridDim = as.integer(gridDim),
        maxdimension = as.integer(maxdimension), decomposition = "5tetrahedra",
		library = library, location = location, printProgress = printProgress)
  } else {
    gridOut <- GridDiag(
	    FUNvalues = FUNvalues, gridDim = as.integer(gridDim),
        maxdimension = as.integer(maxdimension), decomposition = "barycenter",
		library = library, location = location, printProgress = printProgress)
  }

  if (location == TRUE) {
    BirthLocation <- Grid[["grid"]][gridOut[[2]][, 1], ]
    DeathLocation <- Grid[["grid"]][gridOut[[2]][, 2], ]
    if (library == "Dionysus") {
      CycleLocation <- lapply(gridOut[[3]], function(bdy) {
          array(Grid[["grid"]][bdy, ], dim = c(dim(bdy), NCOL(Grid[["grid"]])))})
    }
  }

  Diag <- gridOut[[1]]
  if (NROW(Diag) > 0) {
    Diag[1, 3] <- ifelse(is.null(diagLimit), max(FUNvalues), diagLimit) 
  }
  if (sublevel == FALSE) {
    colnames(Diag) <- c("dimension", "Death", "Birth")
    Diag[, 2:3] <- -Diag[, 3:2]
  } else {
    colnames(Diag) <- c("dimension", "Birth", "Death")
  }

  class(Diag) <- "diagram"
  attributes(Diag)[["maxdimension"]] <- max(Diag[, 1])
  nonInf <- which(Diag[, 2] != Inf & Diag[, 3] != Inf)
  attributes(Diag)[["scale"]] <-
    c(min(Diag[nonInf, 2:3]), max(Diag[nonInf, 2:3]))
  attributes(Diag)[["call"]] <- match.call()
  if (location == FALSE || library == "GUDHI")
  {
    out <- list("diagram" = Diag)
  } else if (library == "PHAT")
  {
    out <- list("diagram" = Diag, "birthLocation" = BirthLocation,
        "deathLocation" = DeathLocation)
  } else
  {
    out <- list("diagram" = Diag, "birthLocation" = BirthLocation,
        "deathLocation" = DeathLocation, "cycleLocation" = CycleLocation)
  }
  return (out)
}
