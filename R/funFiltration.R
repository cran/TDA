funFiltration <- function(FUNvalues, cmplx, sublevel = TRUE) {

  if (!is.numeric(FUNvalues)) {
    stop("FUNvalues should be a vector")
  }

  if (!is.list(cmplx)) {
    stop("cmplx should be a numeric list")
  }

  if (!is.logical(sublevel)) {
    stop("sublevel should be logical")
  }

  if (sublevel == FALSE) {
    FUNvalues <- -FUNvalues
  }

  funOut <-
      FunFiltration(FUNvalues = FUNvalues, cmplx = cmplx)

  if (sublevel == TRUE) {
    out <- list(
        "cmplx" = funOut[[1]], "values" = funOut[[2]], "increasing" = sublevel)
  } else {
    out <- list(
        "cmplx" = funOut[[1]], "values" = -funOut[[2]], "increasing" = sublevel)
  }

  return (out)
}