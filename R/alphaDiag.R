alphaDiag <-
function(X, library = "GUDHI", printProgress = FALSE) {

  # in 32bit architectures alphaDiag doesn't work
  if (.Machine[["sizeof.pointer"]] != 8) {
    Diag <- matrix(0, nrow = 0, ncol = 3)
    class(Diag) <- "diagram"
    attributes(Diag)[["maxdimension"]] <- -Inf
    attributes(Diag)[["scale"]] <- c(Inf, -Inf)
    attributes(Diag)[["call"]] <- match.call()
    out <- list("diagram" = Diag)
    cat("alphaDiag function currently only work on 64-bit R")
    return (out)
  }

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
  if (library == "gudhi" || library == "Gudhi") {
    library <- "GUDHI"
  }
  if (library != "GUDHI") {
    stop("library should be 'GUDHI'")
  }
  if (!is.logical(printProgress)) {
    stop("printProgress should be logical")
  }

  X <- as.matrix(X)

  if (library == "GUDHI") {
    alphaOut <- AlphaDiagGUDHI(X = X, printProgress = printProgress)
  }

  Diag <- alphaOut[[1]]

  colnames(Diag) <- c("dimension", "Birth", "Death")
  class(Diag) <- "diagram"
  attributes(Diag)[["maxdimension"]] <- max(Diag[, 1])
  nonInf <- which(Diag[, 2] != Inf & Diag[, 3] != Inf)
  attributes(Diag)[["scale"]] <-
      c(min(Diag[nonInf, 2:3]), max(Diag[nonInf, 2:3]))
  attributes(Diag)[["call"]] <- match.call()
  out <- list("diagram" = Diag)
  return (out)
}
