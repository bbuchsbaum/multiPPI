#include <RcppArmadillo.h>
using namespace Rcpp;

// [[Rcpp::export]]
arma::mat ridge_precision_cpp(const arma::mat& S, double lambda) {
  const int p = S.n_rows;
  if (S.n_cols != p) stop("S must be square.");
  arma::mat Sr = S + lambda * arma::eye(p, p);
  return arma::inv_sympd(Sr);
}

// [[Rcpp::export]]
arma::cube precision_gate_cpp(const arma::mat& Theta0, const arma::cube& M) {
  const int r = Theta0.n_rows;
  if (Theta0.n_cols != r) stop("Theta0 must be square.");
  if (M.n_rows != r || M.n_cols != r)
    stop("Cube dimensions do not align with Theta0.");
  const int K = M.n_slices;
  arma::cube out(r, r, K);
  for (int k = 0; k < K; ++k) {
    out.slice(k) = -Theta0 * M.slice(k) * Theta0;
  }
  return out;
}
