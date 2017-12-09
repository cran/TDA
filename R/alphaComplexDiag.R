alphaComplexDiag <-
# function(X, maxdimension = NCOL(X) - 1, maxalphasquare, library = "GUDHI",
#     location = FALSE, printProgress = FALSE) {
function(X, maxdimension = NCOL(X) - 1, library = "GUDHI", location = FALSE,
    printProgress = FALSE) {

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (!is.numeric(maxdimension) || length(maxdimension) != 1 ||
      maxdimension < 0) {
    stop("maxdimnsion should be a nonnegative integer")
  }
#  if (!is.numeric(maxalphasquare) || length(maxalphasquare) != 1 || maxalphasquare < 0) {
#    stop("maxalphasquare should be a nonnegative number")
#  }
  if (length(library) == 1) {
    library <- rep(library, 2)
  }
  if (library[1] == "gudhi" || library[1] == "Gudhi") {
    library[1] <- "GUDHI"
  }
  if (library[1] != "GUDHI") {
    stop("library for building a filtration should be 'GUDHI'")
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

  X <- as.matrix(X)
  maxdimension <- min(maxdimension, NCOL(X) - 1)

#  alphaOut <- AlphaComplexDiag(
#      X = X, maxalphasquare = maxalphasquare, libraryDiag = library[2],
#      location = location, printProgress = printProgress)  
  alphaOut <- AlphaComplexDiag(
   X = X, maxdimension = maxdimension, libraryDiag = library[2],
      location = location, printProgress = printProgress)

  if (location == TRUE) {
    BirthLocation <- X[alphaOut[[2]][, 1], ]
    DeathLocation <- X[alphaOut[[2]][, 2], ]
    if (library[2] == "Dionysus") {
      CycleLocation <- lapply(alphaOut[[3]], function(bdy) {
          array(X[bdy, ], dim = c(dim(bdy), NCOL(X)))})
    }
  }  
  
  Diag <- alphaOut[[1]]

  colnames(Diag) <- c("dimension", "Birth", "Death")
  class(Diag) <- "diagram"
  attributes(Diag)[["maxdimension"]] <- max(Diag[, 1])
  nonInf <- which(Diag[, 2] != Inf & Diag[, 3] != Inf)
  attributes(Diag)[["scale"]] <-
      c(min(Diag[nonInf, 2:3]), max(Diag[nonInf, 2:3]))
  attributes(Diag)[["call"]] <- match.call()

  if (location == FALSE || library[2] == "GUDHI") {
    out <- list("diagram" = Diag)
  } else if (library[2] == "PHAT") {
    out <- list("diagram" = Diag, "birthLocation" = BirthLocation,
        "deathLocation" = DeathLocation)
  } else {
    out <- list("diagram" = Diag, "birthLocation" = BirthLocation,
        "deathLocation" = DeathLocation, "cycleLocation" = CycleLocation)
  }

  return (out)
}
