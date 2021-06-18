whiteNoiseSimmap <- function (Tree,
                              ThresholdTime) {
  Bt <- branching.times(Tree)
  N <- Ntip(Tree)
  Nh <- nodeHeights(Tree)
  Nh <- max(Nh) - Nh
  MaxBt <- max(Bt)
  Slice <- MaxBt + MaxBt/1e6 - ThresholdTime
  RescaleValue <- MaxBt - Slice
  TipSliced <- treeSlice(Tree, slice = Slice, trivial = TRUE)
  More2Tips <- unlist(lapply(TipSliced, function(x) Ntip(x) > 1))
  if (any(More2Tips)) {
    TipSliced <- TipSliced[More2Tips]
    Tips <- lapply(TipSliced, function(x) x$tip.label)
    for (z in 1:length(Tips)) {
      Mrca <- getMRCA(Tree, Tips[[z]])
      W <- Tree$edge[, 2] == Mrca
      Tree$edge.length[W] <- Nh[W, 1] -  RescaleValue
      Desc <- getDescendants(Tree, Mrca)
      Desc <- Desc[Desc > N]
      Tree$edge.length[Tree$edge[, 2] %in% Desc] <- 0
      TipIdx <- which(Tree$tip.label %in% Tips[[z]])
      Tree$edge.length[Tree$edge[, 2] %in% TipIdx] <- RescaleValue
    }
  }
  TreeSimmap <- make.era.map(Tree, c(0, Slice - Slice/1e3))
  return(TreeSimmap)
}
