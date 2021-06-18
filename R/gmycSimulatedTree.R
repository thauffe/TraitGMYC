#' @title Simulate tree for GMYC
#'
#' @description This function is a fork of GMYC.simulated.tree from the P2C2M.GMYC package
#' and graft the speciation and coalescent tree and correct the branch lengths.
#'
#' @usage gmycSimulatedTree(SpeciationTree, SpeciesResult)
#'
#' @param SpeciationTree tree
#' @param SpeciesResult a data.frame with an index for species in the first column,
#' the number of individuals per species in the 2nd column
#' and the nucleotide diversity (as proxy for population size) in the 3rd column
#' (NA for species with one individual)
#' @param Scale should tree be scaled to unit deepth? (default TRUE)
#'
#' @author Torsten Hauffe (with code forked from Emanuel M. Fonseca and Drew J. Duckett)
#'
#' @export gmycSimulatedTree
gmycSimulatedTree <- function (SpeciationTree, SpeciesResult, Scale = TRUE) {
  n.tips <- sum(SpeciesResult$N.samples)
  if (length(SpeciesResult[, 1]) != n.tips & length(SpeciesResult[, 1]) != 1) {
    number <- 0
    for (u in 1:nrow(SpeciesResult)) {
      if (SpeciesResult[u, 2] > 1) {
        if (number == 0) {
          number <- u
          sim <- sim.coaltree(nspecies = SpeciesResult[u, 2], theta = SpeciesResult[u, 3])
          sim <- read.tree(text = paste(sim, ";", sep = ""))
          sim$tip.label <- paste(rep(paste("Tip", SpeciesResult[u, 1], ".", sep = ""),
                                     e = SpeciesResult[u, 2]), 1:SpeciesResult[u, 2], sep = "")
          if (any(SpeciesResult[, 1] > 1)) {
            Grafted.tree <- bind.tree(SpeciationTree, sim,
                                      where = which(SpeciationTree$tip.label == paste("Tip", u, ".1", sep = "")))
          }
          Grafted.tree <- force.ultrametric(Grafted.tree, method = "nnls")
          if (Scale) {
            Grafted.tree$edge.length <- Grafted.tree$edge.length/max(nodeHeights(Grafted.tree))
          }
        }
        else {
          sim <- sim.coaltree(nspecies = SpeciesResult[u, 2], theta = SpeciesResult[u, 3])
          sim <- read.tree(text = paste(sim, ";", sep = ""))
          sim$tip.label <- paste(rep(paste("Tip", SpeciesResult[u, 1], ".", sep = ""),
                                     e = SpeciesResult[u, 2]), 1:SpeciesResult[u, 2], sep = "")
          Grafted.tree <- bind.tree(Grafted.tree, sim,
                                    where = which(Grafted.tree$tip.label == paste("Tip", u, ".1", sep = "")))
          Grafted.tree <- force.ultrametric(Grafted.tree, method = "nnls")
          if (Scale) {
            Grafted.tree$edge.length <- Grafted.tree$edge.length/max(nodeHeights(Grafted.tree))
          }
        }
      }
    }
  }
  else if (length(SpeciesResult[, 1]) == n.tips) {
    Grafted.tree <- SpeciationTree
  }
  else {
    sim <- sim.coaltree(nspecies = SpeciesResult[1, 2],
                        theta = SpeciesResult[1, 3])
    sim <- read.tree(text = paste(sim, ";", sep = ""))
    sim$tip.label <- paste(rep(paste("Tip", SpeciesResult[1, 1], ".", sep = ""),
                               e = SpeciesResult[1, 2]), 1:SpeciesResult[1, 2], sep = "")
    Grafted.tree <- sim
  }
  return(Grafted.tree)
}







