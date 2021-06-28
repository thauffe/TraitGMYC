flattenMat <- function (Mat, Prefix) {
  Flat <- c(Mat)
  Names <- colnames(Mat)
  Names <- paste0(rep(Names, each = length(Names)), "_", Names)
  Names <- paste0(Prefix, Names)
  names(Flat) <- Names
  return(Flat)
}

fillFail <- function(TraitMat) {
  Traits <- colnames(TraitMat)
  LenTraits <- ncol(TraitMat)
  Mat <- matrix(NA_integer_, LenTraits, LenTraits)
  colnames(Mat) <- rownames(Mat) <- Traits
  Thetas <- rep(NA_integer_, LenTraits)
  names(Thetas) <- paste0("Theta_", Traits)
  Sigmas <- c(flattenMat(Mat, "SigmaSp_"),
              flattenMat(Mat, "SigmaPop_"))
  Out <- c(Thetas, Sigmas)
  return(Out)
}

fitTraitEvolution <- function (ThresholdTime,
                               tr,
                               trait,
                               meserr,
                               maxit,
                               quiet,
                               TraitModel,
                               Method) {
  if (TraitModel == "BMBM") {
    Bt <- branching.times(tr)
    TreeSimmap <- make.era.map(tr, c(0, max(Bt) - ThresholdTime))
  }
  if (TraitModel == "BMWN") {
    TreeSimmap <- whiteNoiseSimmap(tr, ThresholdTime)
  }
  Param <- NULL
  if (ncol(trait) > 1) {
    Param <- list(constraint = "proportional", decomp = "spherical")
  }
  TraitEvo <- tryCatch(mvBM(tree = TreeSimmap,
                            data = trait, error = meserr,
                            model = "BMM",
                            method = Method,
                            param = Param,
                            control = list(maxit = maxit),
                            diagnostic = !quiet,
                            echo = !quiet),
                       error = function(e) NA)
  if (!is.na(TraitEvo[[1]])) {
    Thetas <- c(TraitEvo$theta)
    names(Thetas) <- paste0("Theta_", colnames(TraitEvo$theta))
    Res <- c(ThresholdTime,
             TraitEvo$LogLik,
             Thetas,
             flattenMat(TraitEvo$sigma[,,1], "SigmaSp_"),
             flattenMat(TraitEvo$sigma[,,2], "SigmaPop_"))
  }
  else {
    Res <- c(NA, NA, fillFail(trait))
  }
  names(Res)[1:2] <- c("ThresholdTime", "likelihood")
  return(Res)
}

runTimeSliceMultivariate <- function (gmycResults,
                                      tr,
                                      trait,
                                      meserr,
                                      TraitModel,
                                      maxit = 20000,
                                      quiet = TRUE,
                                      ncores = 1) {
  N <- length(gmycResults$threshold.time)
  Ntraits <- ncol(trait)
  OutputMat <- as.data.frame(matrix(NA_real_,
                                    nrow = N,
                                    ncol = 2 + Ntraits + 2*Ntraits*Ntraits))
  Method <- "rpf"
  if (ncol(trait) == 1 && sum(is.na(trait[, 1])) > 0) {
    Method <- "inverse"
  }
  TraitEvo <- mvBM(tree = tr,
                   data = trait, error = meserr,
                   model = "BM1",
                   method = Method,
                   control = list(maxit = maxit),
                   diagnostic = !quiet,
                   echo = !quiet)
  Thetas <- c(TraitEvo$theta)
  names(Thetas) <- paste0("Theta_", colnames(TraitEvo$theta))
  FitNoShift <- c(max(-gmycResults$threshold.time),
                  TraitEvo$LogLik,
                  Thetas,
                  flattenMat(TraitEvo$sigma, "SigmaSp_"))
  OutputMat[1, 1:length(FitNoShift)] <- FitNoShift
  registerDoParallel(ncores)
  FitTraits <- foreach(iter = 2:N,
                       .packages = c("mvMORPH"),
                       .combine = rbind,
                       .inorder = TRUE) %dopar% fitTraitEvolution(ThresholdTime = -gmycResults$threshold.time[iter],
                                                                  tr = tr,
                                                                  trait = trait,
                                                                  meserr = meserr,
                                                                  maxit = maxit,
                                                                  quiet = quiet,
                                                                  TraitModel = TraitModel,
                                                                  Method = Method)
  stopImplicitCluster()
  OutputMat[2:nrow(OutputMat), ] <- FitTraits
  colnames(OutputMat) <- colnames(FitTraits)
  # for (i in 1:N) {
  #   ThresholdTime <- -gmycResults$threshold.time[i]
  #   if (i == 1) {
  #     TraitEvo <- mvBM(tree = tr,
  #                      data = trait, error = meserr,
  #                      model = "BM1",
  #                      method = "rpf",
  #                      control = list(maxit = maxit),
  #                      diagnostic = !quiet,
  #                      echo = !quiet)
  #     OutputMat[i, 1:(2+ncol(trait))] <- c(ThresholdTime,
  #                                          TraitEvo$LogLik,
  #                                          diag(TraitEvo$sigma))
  #   }
  #   else {
  #     Bt <- branching.times(tr)
  #     TreeSimmap <- make.era.map(tr, c(0, max(Bt) - ThresholdTime))
  #     TraitEvo <- tryCatch(mvBM(tree = TreeSimmap,
  #                               data = trait, error = meserr,
  #                               model = "BMM",
  #                               method = "rpf",
  #                               param = list(constraint = "proportional", decomp = "eigen+"),
  #                               control = list(maxit = maxit),
  #                               diagnostic = !quiet,
  #                               echo = !quiet),
  #                          error = function(e) NA)
  #     if (!is.na(TraitEvo[[1]])) {
  #       OutputMat[i, ] <- c(ThresholdTime,
  #                           TraitEvo$LogLik,
  #                           .flattenSigma(TraitEvo$sigma[,,1]),
  #                           .flattenSigma(TraitEvo$sigma[,,2]))
  #     }
  #   }
  # }
  return(OutputMat)
}
