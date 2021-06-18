.getEdgePosition <- function(H, thresh) {
  if (H[2] >= thresh) {
    edgePos <- "before"
  }
  else if (H[1] > thresh)  {
    edgePos <- "crosses"
  }
  else {
    edgePos <- "after"
  }
  return(edgePos)
}
