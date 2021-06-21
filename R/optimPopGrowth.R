harmonicMean <- function (x) {
  1/((1/length(x))*sum(1/x))
}
nT <- function (N0, Alpha, TGrowth) {
  N0 * exp(-Alpha * (seq(0, TGrowth, TGrowth/100)))
}
optimFunction <- function (N0, Alpha, TGrowth, PopSize) {
  H <- harmonicMean( nT(N0, Alpha, TGrowth) )
  O <- (H - PopSize)^2
  return(O)
}