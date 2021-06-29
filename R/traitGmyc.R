#' @title Optimizes genetic and phenotypic
#' clusters using a combination of the generalized mixed yule coalescent and
#' models of trait evolution.
#'
#' @description This function fits within- and between-species branching models
#' to reconstructed gene trees, known as the generalized mixed yule coalescent (GMYC)
#' model, and the Brownian motion process to (continous) trait evolution with
#' a change in evolutionary rates.
#'
#' @param tr An ultrametric, dichotomous tree object in ape format.
#' @param interval Upper and lower limit of estimation of scaling parameters, e.g. c(0,10)
#' @param trait A matrix of trait values with rownames corresponding to tips of the tree (NA for missing traits)
#' @param traitmodel Model of trait evolution.
#' Either "BMBM" where a shift in rate of morphological variation is modeled.
#' Or "BMWN" where a shift in rate together with a shift in mode of evolution
#' from Brownian motion towards no phylogenetic signal is modeled (i.e. white noise).
#' @param quiet By default shows no progress on console. Use quiet = TRUE to enable.
#' @param meserr A data.frame (or matrix) of squared standard error for each traits.
#'               Can contain NAs. row.names should match tip labels of the phylogeny.
#' @param ncores Number of cores used for fitting models of trait evolution
#'
#' @return traitgmyc returns an object of class "traitgmyc": a list with the following elements
#' \item{method}{ method used for an analysis}
#' \item{likelihood}{ likelihood values for each gmyc optimization}
#' \item{parameters}{ estimated parameters for each gmyc optimization. (lambda1, lambda2, pp1, pp2)}
#' \item{entity}{ numbers of entities}
#' \item{cluster}{ numbers of clusters}
#' \item{MRCA}{ index of MRCA nodes, i.e. ancestral node of each delimited cluster}
#' \item{threshold.time}{ optimized threshold times}
#' \item{tree}{ the tree}
#' \item{traitsGmyc}{ data.frame with likelihood and rate parameters for Brownian motion process with a shift at gmyc threshold times}
#' \item{trait}{ traits}
#' \item{sum_likelihoods}{ sum of gmyc and Brownian motion likelihood}
#'
#' @author Torsten Hauffe and Robin Ackermann
#'
#' @examples
#' \dontrun{
#' N <- 10
#' SpeciesResult <- data.frame(Species = 1:N,
#'                             Individuals = sample(2:20, N, replace = TRUE))
#' SpeciesResult$Theta <- runif(N, min = 0.01, max = 0.7)
#' Tree <- pbtree(b = 0.2, n = N)
#' Tree <- replaceTiplabel(Tree, Newlabel = "Tip")
#' GmycTree <- gmycSimulatedTree(Tree, SpeciesResult, Scale = FALSE)
#' GmycTreePainted <- paintSpeciesBranches(GmycTree)
#'
#' Ntraits <- 3
#' SigmasSpecies <- simSigma(Ntraits)
#' Cor <- cov2cor(SigmasSpecies)
#' PopSigmaMulti <- 2
#' SigmasPopulations <- simSigma(Ntraits,
#'                               Cor = Cor[lower.tri(Cor)],
#'                               Sigma2 = PopSigmaMulti * sqrt(diag(SigmasSpecies)))
#'
#' Sigmas <- list(Species = SigmasSpecies, Populations = SigmasPopulations)
#' SimTraits <- mvSIM(GmycTreePainted, model = "BMM",
#'                    param = list(ntraits = Ntraits,
#'                                 theta = rep(0, Ntraits),
#'                                 sigma = Sigmas))
#' SimTraits[1,1] <- NA
#'
#' Res <- traitGmyc(tr = GmycTree,
#'                  trait = SimTraits,
#'                  meserr = NULL,
#'                  quiet = TRUE,
#'                  ncores = 1)
#' plot(Res)
#' 
#' SpeciesTree <- pbtree(b = 0.27, n = 10)
#' GeneTree <- simGenealogy(Species = SpeciesTree,
#'                          Scenario = "B",
#'                          Ind = 5,
#'                          PopSize = 100000)
#' SimTraits <- simTraitsIndividuals(SpeciesTree,
#'                                   Ntraits = 4,
#'                                   IndPop = GeneTree$Species,
#'                                   Sigma2 = rep(1, 4))
#' Res <- traitGmyc(tr = GeneTree$Genealogy,
#'                  trait = SimTraits,
#'                  traitmodel = "BMWN",
#'                  meserr = NULL,
#'                  quiet = TRUE,
#'                  ncores = 1)
#' plot(Res, ask = FALSE)
#' }

#' @export traitGmyc
traitGmyc <- function (tr,
                       interval = c(0, 5),
                       trait,
                       meserr = NULL,
                       traitmodel = "BMBM",
                       quiet = TRUE,
                       ncores = 1) {
  if (quiet) {
    Sink <- capture.output({
      gmycResults <- gmyc(tr = tr,
                          interval = interval,
                          quiet = quiet)
    })
  }
  else {
    gmycResults <- gmyc(tr = tr,
                        interval = interval,
                        quiet = quiet)
  }
  trait <- trait[pmatch(tr$tip.label, rownames(trait)), , drop = FALSE]
  meserr <- meserr[pmatch(tr$tip.label, rownames(meserr)), , drop = FALSE]
  if (is.null(colnames(trait))) {
    colnames(trait) <- paste0("Trait", 1:ncol(trait))
  }
  traitResults <- runTimeSliceMultivariate(gmycResults = gmycResults,
                                           tr = tr,
                                           trait = trait,
                                           meserr = meserr,
                                           TraitModel = traitmodel,
                                           quiet = quiet,
                                           ncores = ncores)
  gmycResults[["traitsGmyc"]] <- traitResults
  gmycResults[["trait"]] <- trait
  gmycResults[["meserr"]] <- meserr
  gmycResults[["sum_likelihoods"]] <- gmycResults$likelihood + traitResults$likelihood
  gmycResults[["TraitModel"]] <- traitmodel
  class(gmycResults) <- "traitgmyc"
  return(gmycResults)
  return(gmycResults)
}
