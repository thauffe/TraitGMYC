#' @title Simulate tree for GMYC
#'
#' @description This function is a fork of GMYC.simulated.tree from the P2C2M.GMYC package
#' It simulates a tree under the yule (pure-birth) process with the same number of species
#' that was delimited using the empirical dataset
#'
#' @usage gmycSimulatedTree(SpeciationTree, EmpiricalTree)
#'
#' @param SpeciesResult a data.frame with an index for species in the first column,
#' the number of individuals per species in the 2nd column
#' and the nucleotide diversity (as proxy for population size) in the 3rd column
#' (NA for species with one individual)
#' @param EmpiricalTree Empirical tree used for species delimitation
#' @param Scale should tree be scaled to unit depth? (default TRUE)
#'
#' @author Torsten Hauffe (with code forked from Emanuel M. Fonseca and Drew J. Duckett)
#'
#' @export yuleTree


yuleTree <- function (SpeciesResult, EmpiricalTree, Scale = TRUE)
{
  tips <- EmpiricalTree$tip.label
  TreeDropped <- drop.tip(EmpiricalTree, tips, trim.internal = FALSE)
  D1 <- max(nodeHeights(EmpiricalTree))
  D2 <- max(nodeHeights(TreeDropped))
  Add <- D1 - D2
  TreeDropped  <- extantBranches2Present(tree = TreeDropped , tol = 1e-08, Add = Add)
  YuleFit <- yule(TreeDropped)
  yule_lambda <- rnorm(1, YuleFit$lambda, YuleFit$se)
  if (nrow(SpeciesResult) > 1) {
    yule.tree <- pbtree(b = yule_lambda, n = N)
    yule.tree <- replaceTiplabel(yule.tree, Newlabel = "Tip")
    yule.tree$edge.length[yule.tree$edge[, 2] <= Ntip(yule.tree)] <- yule.tree$edge.length[yule.tree$edge[, 2] <= Ntip(yule.tree)] + 0.05
    if (Scale) {
      yule.tree$edge.length <- yule.tree$edge.length/max(nodeHeights(yule.tree))
    }
    return(yule.tree)
  }
  else {
    return(NULL)
  }
}
