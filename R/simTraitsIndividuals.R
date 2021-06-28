#' @title Simulate traits
#'
#' @description This function simulates traits for a given number of
#' individuals per species and the population size for the species
#'
#' @param Tree The species tree
#' @param Ntraits Number of traits to simulate
#' @param IndPop Matrix or data.frame with the number of individuals
#' per species in the first column and the species' population size in the second column.
#' This element is returned by simGenealogy ('Species').
#' @param PopFluc Optional matrix or data.frame with the parameters for a fluctuating population size
#' (returned as third element form simGenealogy)
#' @param GenLen Optional vector of length one for generation length (the same for all species)
#' @param Sigma2 Optional Brownian motion rate for trait evolution.
#' Vector of length one (the same for all traits and species) or a seperate rate per species.
#' @param Cor Optional correlation strength between traits in the range [-1,1]. 
#' Either a vector of length one (same strength between all traits) or Ntraits*(Ntraits-1)/2.
#' @param Alpha Optional strength of pull towards the morphological optima.
#' Only effective for scenarios A, B and C. 
#' If not specified, optima is (theoretically) reached after half of the tree depth.
#' @parama Scenario Scenario of Fujisawa and Barraclough (2013)
#' * A Null model assuming a neutral coalescent process in a single population
#' * B Diversification (coalescence within a species tree)
#' * D1 Fluctuating population size; bottleneck and then exponential growth
#' * D2 Fluctuating population size; instant growth then shrink
#' * E Diversification with different sized populations
#' * F1 Random sample of individuals per species
#' * F2 Random sample of individuals per species, with sampling probabilities proportional to population sizes
#' * G Geographic structure
#' @param TraitModel Model of trait evolution. Ether Brownian motion (BM)
#' or Orenstein-Uhlenbeck (OU).
#' @md
#' 
#' @author Torsten Hauffe
#'
#' @return A matrix with traits in columns and names of individuals per species as rownames
#'
#' @references Fujisawa, T. and T. Barraclough (2013): Delimiting species using
#' single-locus data and the Generalized Mixed Yule Coalescent approach:
#' a revised method and evaluation on simulated data sets.
#' Systematic Biology 62(5), 707-724.
#'
#' @export simTraitsIndividuals
#'
#'@examples

simTraitsIndividuals <- function(Tree,
                                 Ntraits = 1,
                                 IndPop,
                                 PopFluc = NULL,
                                 GenLength = 1,
                                 Sigma2 = NULL,
                                 Cor = NULL,
                                 Alpha = NULL,
                                 Scenario = "A",
                                 TraitModel = "BM") {
  if (Scenario %in% c("A", "B", "C1")) {
    if ( !all(IndPop[, 2] == IndPop[1, 2]) ) {
      stop("Population sizes should be equal for all species")
    }
    if (Ntraits > 1 && !is.null(Sigma2)) {
      Sigma2 <- sqrt(Sigma2)
    }
    Tol <- 1e-6
    Sigma <- simSigma(Ntraits, Cor = Cor, Sigma2 = Sigma2)
    eS <- eigen(Sigma, symmetric = TRUE)
    ev <- eS$values
    if ( !all( ev >= -Tol * abs(ev[1L]) ) ) {
      Sigma <- as.matrix(nearPD(Sigma)$mat)
    }
    
    if (is.null(Alpha) && TraitModel == "OU") {
      MaxBt <- max(branching.times(Tree))
      # Phylogenetic half-life
      # t1/2 equal to the height of the phylogeny is a moderate value (Cooper et al. 2016)
      A <- log(2) / MaxBt
      A <- rep(Alpha, Ntraits)
      Alpha <- matrix(NA_real_, Ntraits, Ntraits)
      diag(Alpha) <- A
    }
    TraitsSim <- mvSIM(Tree, nsim = 1,
                       model = paste0(TraitModel, 1),
                       param = list(theta = rep(0, Ntraits),
                                    sigma = Sigma,
                                    alpha = Alpha))
    
  }
  # else {
  #   simTraitsMicroEvo()
  # }
  IndTraitSpecies <- indTraitFromSpecies(TraitMeans = TraitsSim,
                                         IndPop = IndPop,
                                         Sigma = Sigma,
                                         GenLength = GenLength)
  return(IndTraitSpecies)
}