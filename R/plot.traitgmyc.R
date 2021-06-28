#' @title Plot the result of traitgmyc
#'
#' @description This function... .
#'
#' @param x An object of class traitgmyc.
#' @param ask logical. If TRUE (and the R session is interactive) the user is asked for input, before a new figure is drawn.
#' @param ... Other parameters to be passed through to plotting functions.
#'
#' @return A plot with three panels.
#'
#' @author Torsten Hauffe and Robin Ackermann
#'
#' @examples
#' 1 + 2
#' # plot()

plot.traitgmyc <- function (x, ask = TRUE, ...) {
  devAskNewPage(ask)
  whiLike <- which.max(x$sum_likelihoods)
  thresh <- -x$threshold.time[whiLike]
  likDiffSmaller2 <- abs(x$sum_likelihoods - max(x$sum_likelihoods, na.rm = TRUE)) < 2
  threshRange <- range(-x$threshold.time[likDiffSmaller2], na.rm = TRUE)
  tree <- x$tree
  Bt <- sort(branching.times(tree), decreasing = TRUE)
  if (x$TraitModel == "BMBM") {
    TreeSimmap <- make.era.map(tree, c(0, max(Bt) - thresh))
  }
  if (x$TraitModel == "BMWN") {
    TreeSimmap <- whiteNoiseSimmap(tree, thresh)
  }
  trait <- x$trait
  Param <- NULL
  if (ncol(trait) > 1) {
    Param <- list(constraint = "proportional", decomp = "spherical")
  }
  Method <- "rpf"
  if (ncol(trait) == 1 && sum(is.na(trait[, 1])) > 0) {
    Method <- "inverse"
  }
  TraitModel <- "BMM"
  if (whiLike == 1) {
    TraitModel <- "BM1"
    Param <- NULL
    TreeSimmap <- tree
  }
  TraitEvo <- mvBM(tree = TreeSimmap,
                   data = trait, error = x$meserr,
                   model = TraitModel,
                   method = Method,
                   param = Param,
                   control = list(maxit = 10000000),
                   diagnostic = FALSE,
                   echo = FALSE)
  if (any(is.na(x$trait))) {
    Imput <- estim(TreeSimmap, data = x$trait, object = TraitEvo,
                   error = x$meserr, asr = FALSE)
    trait <- Imput$estimates
  }
  Asr <- estim(TreeSimmap, data = trait, object = TraitEvo,
               error = x$meserr, asr = TRUE)
  TraitAsr <- rbind(trait, Asr$estimates)
  if (ncol(TraitAsr) > 1) {
    TraitPCA <- prcomp(TraitAsr, scale = FALSE)
    TraitAsr <- TraitPCA$x[, 1]
  }
  Y <- matrix(c(TraitAsr)[tree$edge], nrow(tree$edge), 2)
  Lineages <- 1 + (1:length(Bt))
  plot(Bt, Lineages, type = "n",
       xlim = c(max(-x$threshold.time), 0),
       log = "y", xlab = "Time", ylab = "N")
  rect(xleft = threshRange[2], ytop = length(Bt)*1.2,
       xright = threshRange[1], ybottom = 0.1,
       col = adjustcolor("red", 0.3), border = NA)
  abline(v = thresh, lty = 2, col = "red")
  lines(Bt, Lineages, type = "s")
  plot(-x$threshold.time, x$sum_likelihoods, type = "n",
       xlim = c(max(-x$threshold.time), 0),
       xlab = "Time", ylab = "Likelihood")
  rect(xleft = threshRange[2], ytop = par("usr")[4],
       xright = threshRange[1], ybottom = par("usr")[3],
       col = adjustcolor("red", 0.3), border = NA)
  abline(v = thresh, lty = 2, col = "red")
  lines(-x$threshold.time, x$sum_likelihoods)
  # Plot phenogram with time from the past to the present
  H <- nodeHeights(tree)
  H <- max(H) - H
  edgePos <- apply(H, 1, function(x) .getEdgePosition(x, thresh))
  plot.new()
  plot.window(xlim = c(max(-x$threshold.time), 0), ylim = range(Y))
  axis(1)
  axis(2)
  for(i in 1:nrow(tree$edge)){
    if (edgePos[i] == "before") {
      lines(x = H[i,], y = Y[i,], col = "black", ...)
    }
    else if (edgePos[i] == "crosses")  {
      Y2 <- approx(x = H[i, ], y = Y[i,], xout = thresh)$y
      lines(x = c(H[i, 1], thresh), y = c(Y[i, 1], Y2), col = "black", ...)
      lines(x = c(thresh, H[i, 2]), y = c(Y2, Y[i, 2]), col = "red", ...)
    }
    else {
      lines(x = H[i,], y = Y[i,], col = "red", ...)
    }
  }
  title(xlab = "Time", ylab = "Phenotype")
  layout(matrix(1))
  devAskNewPage(ask = FALSE)
}
