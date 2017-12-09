ripsFiltration <- function(
    X, maxdimension, maxscale, dist = "euclidean", library = "GUDHI",
	printProgress = FALSE) {

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.numeric(maxdimension) ||
      length(maxdimension) != 1 || maxdimension < 0) {
    stop("maxdimension should be a nonnegative integer")
  }
  if (!is.numeric(maxscale) || length(maxscale) != 1) {
    stop("maxscale should be a number")
  }
  if (dist != "euclidean" && dist != "arbitrary") {
    stop ("dist should be either 'euclidean' or 'arbitrary'")
  }
  if (library == "gudhi" || library == "Gudhi") {
    library <- "GUDHI"
  }
  if (library == "dionysus" || library == "DIONYSUS") {
    library <- "Dionysus"
  }
  if (library != "GUDHI" && library != "Dionysus") {
    stop("library should be a string: either 'GUDHI' or 'Dionysus'")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }

  if (dist == "arbitrary" && library == "GUDHI") {
    library <- "Dionysus"
  }

  X <- as.matrix(X)

  max_num_pairs <- 5000  # to be added as an option

  # in 32bit architectures Dionysus L2 doesn't work
  #ripsOut <- RipsDiag(X = X, maxdimension = maxdimension,
  #    maxscale = maxscale, dist = dist, library = library,
  #    location = location, printProgress = printProgress)
  if (dist == "euclidean" && library != "GUDHI" &&
      .Machine[["sizeof.pointer"]] != 8) {
    ripsOut <- RipsFiltration(
	    X = as.matrix(dist(X)), maxdimension = maxdimension,
		maxscale = maxscale, dist = "arbitrary", library = library,
        printProgress = printProgress)
  } else {
    ripsOut <- RipsFiltration(
	    X = X, maxdimension = maxdimension, maxscale = maxscale, dist = dist,
        library = library, printProgress = printProgress)
  }
  
  out <- list(
      "cmplx" = ripsOut[[1]], "values" = ripsOut[[2]], "increasing" = TRUE)

  return (out)
}
