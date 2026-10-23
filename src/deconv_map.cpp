// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <algorithm>
#include <cmath>
#include <vector>

// [[Rcpp::export]]
arma::mat mppi_deconv_map(const arma::mat& Y, const arma::vec& h, const double lambda) {
  if (Y.n_rows == 0 || h.n_elem == 0) {
    Rcpp::stop("MAP deconvolution requires at least one observation and one HRF sample.");
  }
  if (!Y.is_finite() || !h.is_finite() || !std::isfinite(lambda) || lambda < 0) {
    Rcpp::stop("MAP deconvolution requires finite data and HRF and a finite non-negative lambda.");
  }
  const arma::uword n = Y.n_rows;
  const arma::uword kernel_length = std::min(h.n_elem, n);
  // H is triangular with diagonal h[0]. For lambda > 0 the only null
  // vector of D is constant, and H identifies it iff the effective HRF
  // contains a nonzero sample. Handle these exact singular cases explicitly.
  bool identifies_constant = false;
  for (arma::uword k = 0; k < kernel_length; ++k) identifies_constant |= (h[k] != 0.0);
  if ((lambda == 0.0 && h[0] == 0.0) || !identifies_constant) {
    Rcpp::stop("MAP deconvolution system is singular or not positive definite.");
  }
  const arma::uword band = std::min(n - 1, std::max(arma::uword(1), kernel_length - 1));
  const arma::uword width = band + 1;
  // Row i stores the lower-triangular entries (i,i), (i,i-1), ... .
  std::vector<double> lower(static_cast<std::size_t>(n) * width, 0.0);
  auto at = [&](arma::uword i, arma::uword j) -> double& {
    return lower[static_cast<std::size_t>(i) * width + (i - j)];
  };

  // Construct the lower band of H' H + lambda D' D, where H is causal
  // convolution and D contains adjacent first differences. A one-row
  // series has an empty D and therefore no penalty.
  for (arma::uword i = 0; i < n; ++i) {
    const arma::uword first = (i > band) ? i - band : 0;
    for (arma::uword j = first; j <= i; ++j) {
      const arma::uword delta = i - j;
      double value = 0.0;
      if (delta < kernel_length) {
        const arma::uword count = std::min(n - i, kernel_length - delta);
        for (arma::uword k = 0; k < count; ++k) value += h[k] * h[k + delta];
      }
      if (n > 1 && delta == 0) value += lambda * ((i == 0 || i == n - 1) ? 1.0 : 2.0);
      if (delta == 1) value -= lambda;
      at(i, j) = value;
    }
  }

  // In-place banded Cholesky. The factor and the system have the same
  // bandwidth; failure is reported without adding a numerical ridge.
  for (arma::uword i = 0; i < n; ++i) {
    const arma::uword first = (i > band) ? i - band : 0;
    for (arma::uword j = first; j <= i; ++j) {
      double value = at(i, j);
      const arma::uword start = std::max(first, (j > band) ? j - band : arma::uword(0));
      for (arma::uword k = start; k < j; ++k) value -= at(i, k) * at(j, k);
      if (i == j) {
        if (!(value > 0.0) || !std::isfinite(value)) {
          Rcpp::stop("MAP deconvolution system is singular or not positive definite.");
        }
        at(i, j) = std::sqrt(value);
      } else {
        at(i, j) = value / at(j, j);
      }
    }
  }

  arma::mat result(n, Y.n_cols, arma::fill::zeros);
  for (arma::uword col = 0; col < Y.n_cols; ++col) {
    // Forward substitution starts with H'Y; each target shares the factor.
    for (arma::uword i = 0; i < n; ++i) {
      double value = 0.0;
      const arma::uword count = std::min(n - i, kernel_length);
      for (arma::uword k = 0; k < count; ++k) value += h[k] * Y(i + k, col);
      const arma::uword first = (i > band) ? i - band : 0;
      for (arma::uword j = first; j < i; ++j) value -= at(i, j) * result(j, col);
      result(i, col) = value / at(i, i);
    }
    for (arma::uword offset = 0; offset < n; ++offset) {
      const arma::uword i = n - 1 - offset;
      double value = result(i, col);
      const arma::uword last = std::min(n - 1, i + band);
      for (arma::uword j = i + 1; j <= last; ++j) value -= at(j, i) * result(j, col);
      result(i, col) = value / at(i, i);
    }
  }
  return result;
}
