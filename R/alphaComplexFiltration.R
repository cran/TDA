alphaComplexFiltration <- function(
#   X, maxalphasquare, library = "GUDHI", printProgress = FALSE) {
    X, library = "GUDHI", printProgress = FALSE) {

  if (!is.numeric(X) && !is.data.frame(X)) {
    stop("X should be a matrix of coordinates")
  }
#  if (!is.numeric(maxalphasquare) || length(maxalphasquare) != 1 || maxalphasquare < 0) {
#    stop("maxalphasquare should be a nonnegative number")
#  }
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
#   alphaOut <- AlphaComplexFiltration(
#       X = X, maxalphasquare = maxalphasquare, printProgress = printProgress)
    alphaOut <- AlphaComplexFiltration(X = X, printProgress = printProgress)
  }

  out <- list(
      "cmplx" = alphaOut[[1]], "values" = alphaOut[[2]], "increasing" = TRUE)

  return (out)
}
