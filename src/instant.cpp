#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

static inline double fetch_with_fill(const NumericVector &vec, int idx, const std::string &fill) {
  const int n = vec.size();
  if (idx >= 0 && idx < n) return vec[idx];
  if (fill == "edge") {
    if (idx < 0) return vec[0];
    return vec[n - 1];
  }
  if (fill == "na") return NA_REAL;
  return 0.0;
}

// [[Rcpp::export]]
NumericVector inst_corr_rcpp(NumericVector x, NumericVector y,
                             double tau_half,
                             int offset = 0, int warmup = -1,
                             std::string fill = "zero") {
  const int n = x.size();
  if (y.size() != n) stop("Length mismatch between x and y.");
  if (!R_finite(tau_half) || tau_half <= 0) stop("tau_half must be positive.");
  const double eta = std::log(2.0) / tau_half;
  const double a = std::exp(-eta);
  const double one_minus_a = 1.0 - a;

  double Sxy = 0.0, Sxx = 0.0, Syy = 0.0;
  NumericVector out(n, NA_REAL);

  for (int t = 0; t < n; ++t) {
    double xv = x[t];
    double yv = fetch_with_fill(y, t - offset, fill);
    if (!R_finite(xv) || !R_finite(yv)) {
      Sxy *= a;
      Sxx *= a;
      Syy *= a;
      continue;
    }
    Sxy = a * Sxy + one_minus_a * (xv * yv);
    Sxx = a * Sxx + one_minus_a * (xv * xv);
    Syy = a * Syy + one_minus_a * (yv * yv);
    if (warmup >= 0 && t < warmup) {
      out[t] = NA_REAL;
    } else {
      out[t] = (Sxx > 0 && Syy > 0) ? Sxy / (std::sqrt(Sxx) * std::sqrt(Syy)) : NA_REAL;
    }
  }

  return out;
}
