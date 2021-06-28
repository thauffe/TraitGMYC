simTraitsMicroEvo <- function (Tree,
                               PopSize = 100,
                               Ntraits = 1,
                               GenLength = 1,
                               Sigma2 = NULL,
                               Cor = NULL) {
  Sigma <- simSigma(Ntraits, Cor = Cor, Sigma2 = Sigma2)
  if (Ntraits > 1) {
    Tol <- 1e-6
    eS <- eigen(Sigma, symmetric = TRUE)
    ev <- eS$values
    if ( !all( ev >= -Tol * abs(ev[1L]) ) ) {
      Sigma <- as.matrix(nearPD(Sigma)$mat)
    }
  }
  MaxBt <- max(branching.times(Tree))
  Sigma <- (Sigma * PopSize) / (1e6/GenLength) # Gives standard deviation

  Tree <- reorder(Tree)
  Tree$edge.length <- Tree$edge.length * (1e6/GenLength)
  NodeTips <- unique(c(Tree$edge))
  NodeTipsMeans <- matrix(NA_real_, ncol = Ntraits, nrow = length(NodeTips))
  rownames(NodeTipsMeans) <- NodeTips
  NodeTipsMeans[1, ] <- rep(0, Ntraits)
  for (i in 1:nrow(Tree$edge)) {
    Time <- round(Tree$edge.length[i])
    FromTo <- Tree$edge[i, ]
    Mu <- NodeTipsMeans[rownames(NodeTipsMeans) == FromTo[1], , drop = FALSE]
    Rcpp <- TRUE
    if (!Rcpp) {
      for (tt in 1:Time) {
        # Traits <- as.matrix(rnorm(n = PopSize, mean = Mu, sd = sqrt(Sigma)), ncol = 1)
        Traits <- mvrnorm(n = PopSize, mu = Mu, Sigma = Sigma)
        Mu <- colMeans(Traits)
      }
    }
    else {
      Mu <- getMvrnormMeansAfterTime(Time, PopSize, Mu, Sigma)
      Mu <- c(Mu)
    }
    NodeTipsMeans[rownames(NodeTipsMeans) == FromTo[2], ] <- Mu
  }
  W <- match(1:Ntip(Tree), NodeTips)
  rownames(NodeTipsMeans)[W] <- Tree$tip.label
  NodeTipsMeans <- NodeTipsMeans[-c(1:(Ntip(Tree)-1)), , drop = FALSE]
  return(NodeTipsMeans)
}
