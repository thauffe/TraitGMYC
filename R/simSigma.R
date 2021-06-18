#' @title Simulate variance-covariance matrix for morphological evolution
#'
#' @description This function generates a variance-covariance matrix for
#' morphological evolution, defining the rate at which traits evolve and their
#' evolutionary correlation.
#'
#' @usage simSigma(Ntraits, Cor = NULL, Sigma2 = NULL)
#'
#' @param Ntraits number of traits
#' @param Cor correlation between traits. Default between -1 and 1.
#' Optional, can be fixed to be equal between all traits by giving one value or
#' Ntraits*(Ntraits-1)/2
#' @param Sigma2 Brownian motion rate. Default between 1e-4 and 0.2
#' Optional, can be fixed to be equal for all traits or one value per trait
#'
#' @author Torsten Hauffe
#'
#' @return matrix Ntrait x Ntrait for simulating trait evolution
#'
#' @export simSigma

simSigma <- function(Ntraits, Cor = NULL, Sigma2 = NULL) {
  if (!is.null(Sigma2)) {
    if ( length(Sigma2) != Ntraits && length(Sigma2) != 1) {
      stop("Sigma2 should be of length 1 or Ntraits")
    }
    else {
      if (length(Sigma2) == 1) {
        Sigma2 <- rep(Sigma2, Ntraits)
      }
    }
  }
  else {
    Sigma2 <- runif(Ntraits, min = 1e-4, max = 0.2)
  }
  Cov <- matrix(1, ncol = Ntraits, nrow = Ntraits)
  Q <- Ntraits*(Ntraits-1)/2
  if (!is.null(Cor)) {
    if (length(Cor) != Q && length(Cor) != 1) {
      stop("Correlation among traits should be of length 1 or Ntraits*(Ntraits-1)/2")
    }
    SimCov <- Cor
  }
  else {
    SimCov <- runif(Q, min = -1, max = 1) # Trait correlation
  }
  Cov[lower.tri(Cov, diag = FALSE)] <- SimCov
  Cov <- t(Cov)
  Cov[lower.tri(Cov, diag = FALSE)] <- SimCov
  Sigmas <- diag(Sigma2)  %*% Cov  %*% diag(Sigma2) # Correlation to covariance
  return(Sigmas)
}

