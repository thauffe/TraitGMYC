// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppArmadillo.h which pulls Rcpp.h in for us
#include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// https://gallery.rcpp.org/articles/simulate-multivariate-normal/
// mu is a vector
// // [[Rcpp::export]]
// arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
//   int ncols = sigma.n_cols;
//   arma::mat Y = arma::randn(n, ncols);
//   return arma::repmat(mu, 1, n).t() + Y * sigma;
// }

// mu is a 1-row matrix
// [[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::mat mu, arma::mat sigma) {
  std::cout << mu;
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


// [[Rcpp::export]]
arma::rowvec colMeansArma(Rcpp::NumericMatrix x){
  arma::mat X = arma::mat(x.begin(), x.nrow(), x.ncol(), false);
  return arma::mean(X, 0);
}

// [[Rcpp::export]]
arma::rowvec meansAfterMvrnormDraw(int n, arma::mat mu, arma::mat sigma){
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  arma::mat S = arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
  return arma::mean(S, 0);
}

// [[Rcpp::export]]
arma::mat getMvrnormMeansAfterTime(int tt, int n, arma::mat mu, arma::mat sigma){
  int ncols = sigma.n_cols;
  arma::mat mu2(ncols, tt, arma::fill::zeros);
  mu2.col(0) = mu.t();
  for (int i = 1; i < tt; i++){
    arma::mat Y = arma::randn(n, ncols);
    arma::mat mu3 = mu2.col(i - 1);
    arma::mat S = arma::repmat(mu3, 1, n).t() + Y * arma::chol(sigma);
    arma::mat mu4 = arma::mean(S, 0);
    mu2.col(i) = mu4.t();
  }
  return mu2.col(tt - 1);
}

// [[Rcpp::export]]
arma::mat notWorking(int tt, int n, arma::mat mu, arma::mat sigma){
  int ncols = sigma.n_cols;
  // I need to initialize all elements before the loop?!
  arma::mat Y = arma::randn(n, ncols);
  arma::mat S = arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
  arma::mat mu2 = arma::mean(S, 0);
  std::cout << mu2; // Comment
  for (int i = 1; i < tt; i++){
    arma::mat Y = arma::randn(n, ncols);
    std::cout << mu2;
    arma::mat S = arma::repmat(mu2, 1, n).t() + Y * arma::chol(sigma);
    arma::mat mu2 = arma::mean(S, 0);
    std::cout << mu2;
  }
  return mu2;
}

// Nearest positive definite matrix in C++
// https://stackoverflow.com/questions/51490499/results-for-calculating-nearest-positive-definite-matrix-are-different-in-r-func
