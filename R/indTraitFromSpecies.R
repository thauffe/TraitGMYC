indTraitFromSpecies <- function(TraitMeans,
                                IndPop,
                                Sigma,
                                GenLength = 1) {
  IndTraitSpecies <- vector(mode = "list", nrow(TraitMeans))
  for (i in 1:nrow(TraitMeans)) {
    SigmaSpecies <- (Sigma * IndPop[i, 2]) / (1e6/GenLength)
    IndTraitSpecies[[i]] <- mvrnorm(n = IndPop[i, 1],
                                    mu = TraitMeans[i, ],
                                    Sigma = SigmaSpecies)
  }
  IndTraitSpecies <- do.call("rbind", IndTraitSpecies)
  rownames(IndTraitSpecies) <- c(sapply(1:nrow(IndPop),
                                        function(x)
                                          paste0(rownames(IndPop)[x],
                                                 ".",
                                                 1:IndPop[x, 1])))
  return(IndTraitSpecies)
}
