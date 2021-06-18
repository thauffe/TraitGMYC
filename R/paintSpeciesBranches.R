paintSpeciesBranches <- function(GmycTree) {
  Taxa <- unlist(lapply(strsplit(GmycTree$tip.label, ".", fixed = TRUE), function(x) x[1]))
  UniqueTaxa <- unique(Taxa)
  for (i in 1:length(UniqueTaxa)) {
    W <- which(Taxa == UniqueTaxa[i])
    if (length(W) > 1) {
      AncNode <- getMRCA(phy = GmycTree, tip = GmycTree$tip.label[W])
      GmycTree <- paintSubTree(tree = GmycTree, node = AncNode,
                               state = "Population", anc.state = "Species")
    }
  }
  return(GmycTree)
}
