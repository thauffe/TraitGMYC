replaceTiplabel <- function (Tree, Newlabel = "Tip") {
  Tree$tip.label <- paste(rep(paste(Newlabel, 1:Ntip(Tree), ".1", sep = "")))
  Tree$tip.label <- sample(Tree$tip.label)
  return(Tree)
}
