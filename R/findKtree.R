findKtree <-
function(bb, parent, children, compBranch, n) {
  kappaTop <- numeric(bb)
  kappaBottom <- numeric(bb)
  for (i in seq_len(bb)) {
    if (is.na(parent[i])) {
      kappaBottom[i] <- NA
      kappaTop[i] <- NA
    } else if (parent[i] == 0) { 
        kappaBottom[i] <- 0 
        if (i <= length(children)) {
          Kchildren <- children[[i]]
          if (!is.null(Kchildren) && !is.na(Kchildren)) {
          kappaTop[i] <-
              (length(compBranch[[i]]) - length(unlist(compBranch[Kchildren]))) / n
          } else {
            kappaTop[i] <- length(compBranch[[i]]) / n
          }
        } else {
          kappaTop[i] <- length(compBranch[[i]]) / n
        }
      } else {
        kappaBottom[i] <- kappaTop[parent[i]]
        if (i <= length(children)) {
          Kchildren <- children[[i]]
          if (!is.null(Kchildren) && length(Kchildren) != 0 && !is.na(Kchildren)) {
            kappaTop[i] <- kappaBottom[i] + (length(compBranch[[i]]) -
                length(unlist(compBranch[Kchildren]))) / n
          } else {
            kappaTop[i] <- kappaBottom[i] + length(compBranch[[i]]) / n
          }
        } else {
          kappaTop[i] <- kappaBottom[i] + length(compBranch[[i]]) / n
        }
      }
    }

  return(list("kappaTop" = kappaTop, "kappaBottom" = kappaBottom))
}
