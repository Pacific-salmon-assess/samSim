// [[Rcpp::depends(RcppArmadillo)]]

#include <RcppArmadillo.h>
using namespace Rcpp;

//' This is a function to rapidly fit linear models and return estimated
//' coefficients.
//'
//' @title Estimate coefficients for linear model.
//' @param y A numeric vector of response variables
//' @param X A numeric matrix of predictor variables
//' @return A numeric list with estimated coefficients and standard deviations.
//' @export
// [[Rcpp::export]]
List fastLm(const arma::vec & y, const arma::mat & X) {

  int n = X.n_rows, k = X.n_cols;

  arma::colvec coef = arma::solve(X, y);
  arma::colvec resid = y - X*coef;

  double sig2 = arma::as_scalar(arma::trans(resid)*resid/(n-k));
  arma::colvec stderrest =
    arma::sqrt(sig2 * arma::diagvec( arma::inv(arma::trans(X)*X)) );

  return List::create(Named("coefficients") = coef,
                      Named("stderr")       = stderrest);
}
