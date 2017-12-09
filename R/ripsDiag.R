ripsDiag <-
function(X, maxdimension, maxscale, dist = "euclidean", library = "GUDHI",
         location = FALSE, printProgress = FALSE) {

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
  if (length(library) == 1) {
    library <- rep(library, 2)
  }
  if (library[1] == "gudhi" || library[1] == "Gudhi") {
    library[1] <- "GUDHI"
  }
  if (library[1] == "dionysus" || library[1] == "DIONYSUS") {
    library[1] <- "Dionysus"
  }
  if (library[1] != "GUDHI" && library[1] != "Dionysus") {
    stop("library for building a filtration should be a string: either 'GUDHI' or 'Dionysus'")
  }
  if (library[2] == "gudhi" || library[2] == "Gudhi") {
    library[2] <- "GUDHI"
  }
  if (library[2] == "dionysus" || library[2] == "DIONYSUS") {
    library[2] <- "Dionysus"
  }
  if (library[2] == "phat" || library[2] == "Phat") {
    library[2] <- "PHAT"
  }
  if (library[2] != "GUDHI" && library[2] != "Dionysus" && library[2] != "PHAT") {
    stop("library for computing persistence diagram should be a string: either 'GUDHI', 'Dionysus', or 'PHAT'")
  }
  if (!is.logical(location)) {
    stop("location should be logical")
  }  
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }

  if (dist == "arbitrary" && library[1] == "GUDHI") {
    library[1] <- "Dionysus"
  }

  X <- as.matrix(X)

  max_num_pairs <- 5000  # to be added as an option

  # in 32bit architectures Dionysus/PHAT L2 doesn't work
  #ripsOut <- RipsDiag(X = X, maxdimension = maxdimension,
  #    maxscale = maxscale, dist = dist, library = library,
  #    location = location, printProgress = printProgress)
  if (dist == "euclidean" && library[1] == "Dionysus" &&
      .Machine[["sizeof.pointer"]] != 8) {
    ripsOut <- RipsDiag(X = as.matrix(dist(X)), maxdimension = maxdimension,
        maxscale = maxscale, dist = "arbitrary",
        libraryFiltration = library[1], libraryDiag = library[2],
        location = location, printProgress = printProgress)
  } else {
    ripsOut <- RipsDiag(X = X, maxdimension = maxdimension,
        maxscale = maxscale, dist = dist, libraryFiltration = library[1],
        libraryDiag = library[2], location = location,
        printProgress = printProgress)
  }

  if (location == TRUE) {
    if (dist == "euclidean") {
      BirthLocation <- X[ripsOut[[2]][, 1], ]
      DeathLocation <- X[ripsOut[[2]][, 2], ]
      if (library[2] == "Dionysus") {
        CycleLocation <- lapply(ripsOut[[3]], function(bdy) {
            array(X[bdy, ], dim = c(dim(bdy), NCOL(X)))})
      }
    } else {
      BirthLocation <- ripsOut[[2]][, 1]
      DeathLocation <- ripsOut[[2]][, 2]
      if (library[2] == "Dionysus")
      {
        CycleLocation <- ripsOut[[3]]
      }
    }
  }

  Diag <- ripsOut[[1]]
  if (NROW(Diag) > 0) {
    ## change Inf values to maxscale
    Diag[which(Diag[, 3] == Inf), 3] <- maxscale  
  }

  colnames(Diag) <- c("dimension", "Birth", "Death")
  class(Diag) <- "diagram"
  attributes(Diag)[["maxdimension"]] <- max(Diag[, 1])
  attributes(Diag)[["scale"]] <- c(0, maxscale)
  attributes(Diag)[["call"]] <- match.call()
  if (location == FALSE || library[2] == "GUDHI")
  {
    out <- list("diagram" = Diag)
  } else if (library[2] == "PHAT")
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
