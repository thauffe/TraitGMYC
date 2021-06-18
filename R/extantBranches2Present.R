extantBranches2Present <- function (tree, tol = 1e-08, Add = 0) {
  Nt <- Ntip(tree)
  H <- nodeHeights(tree)
  tl <- max(H)
  x <- which(H[, 2] <= (tl - tol) & tree$edge[, 2] <= Nt)
  tree$edge.length[x] <- tl - H[x, 1]
  TipEl <- tree$edge[, 2] <= Nt
  tree$edge.length[TipEl] <- tree$edge.length[TipEl] + Add
  return(tree)
}
