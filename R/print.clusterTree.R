print.clusterTree <-
function(x, ...) {
  class(x) <- "list"
  print(x)  
}