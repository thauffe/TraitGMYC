.scaleEdges <- function(tree, H, thresh, edgePos, scaleFacEdg) {
  rescaledEdge <- rep(NA_real_, length(tree$edge.length))
  rescaledEdge[edgePos == "after"] <- tree$edge.length[edgePos == "after"] * scaleFacEdg[2]
  rescaledEdge[edgePos == "before"] <- tree$edge.length[edgePos == "before"] * scaleFacEdg[1]
  rescaledEdge[edgePos == "crosses"] <- (thresh - H[edgePos == "crosses", 2]) * scaleFacEdg[2] +
    (H[edgePos == "crosses", 1] - thresh) * scaleFacEdg[1]
  treeRescaled <- tree
  treeRescaled$edge.length <- rescaledEdge
  return(treeRescaled)
}
