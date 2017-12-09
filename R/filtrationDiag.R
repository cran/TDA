filtrationDiag <- function(
    filtration, maxdimension, library = "GUDHI", location = FALSE,
    printProgress = FALSE, diagLimit = NULL) {

  if (length(filtration[["cmplx"]]) != length(filtration[["values"]])) {
    stop("The length of the simplicial complx should equals the length of the filtration values")
  }
  if (!is.numeric(maxdimension) ||
      length(maxdimension) != 1 || maxdimension < 0) {
    stop("maxdimension should be a nonnegative integer")
  }
  if (library == "gudhi" || library == "Gudhi") {
    library <- "GUDHI"
  }
  if (library == "dionysus" || library == "DIONYSUS") {
    library <- "Dionysus"
  }
  if (library != "GUDHI" && library != "Dionysus") {
    stop("library for computing persistence diagram should be a string: either 'GUDHI' or 'Dionysus'")
  }
  if (!is.logical(location)) {
    stop("location should be logical")
  }  
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }  

  if (filtration[["increasing"]] == FALSE) {
    filtration[["values"]] <- -filtration[["values"]]
  }

  filtrationOut <- FiltrationDiag(
      filtration = filtration, maxdimension = maxdimension, library = library,
	  location = location, printProgress = printProgress)

   if (location == TRUE) {
    BirthLocation <- filtrationOut[[2]][, 1]
    DeathLocation <- filtrationOut[[2]][, 2]
    if (library == "Dionysus") {
      CycleLocation <- filtrationOut[[3]]
    }
  }

  Diag <- filtrationOut[[1]]
  if (NROW(Diag) > 0) {
    Diag[Diag == Inf] <- ifelse(is.null(diagLimit), Inf, diagLimit) 
  }

  if (filtration[["increasing"]] == FALSE) {
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
  if (location == FALSE || library == "GUDHI") {
    out <- list("diagram" = Diag)
  } else if (library == "PHAT") {
    out <- list(
	    "diagram" = Diag, "birthLocation" = BirthLocation,
		"deathLocation" = DeathLocation)
  } else {
    out <- list(
	    "diagram" = Diag, "birthLocation" = BirthLocation,
		"deathLocation" = DeathLocation, "cycleLocation" = CycleLocation)
  }
  return (out)
}