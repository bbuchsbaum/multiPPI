#include <Rcpp.h>
#include <algorithm>
#include <cmath>

using namespace Rcpp;

// [[Rcpp::export]]
List ewm_gated_cross_rcpp(NumericMatrix U,
                          NumericMatrix PK,
                          double eta,
                          IntegerVector offsets,
                          bool normalized = true,
                          IntegerVector blocklens = IntegerVector()) {
  const int T = U.nrow();
  const int r = U.ncol();
  const int K = PK.ncol();
  const int L = offsets.size();

  if (T < 2) stop("Number of time points must be >= 2.");
  if (eta <= 0.0) stop("eta must be positive.");
  if (K == 0 || L == 0) return List::create();

  const double a = std::exp(-eta);
  const double one_minus_a = 1.0 - a;

  std::vector<int> run_id;
  const bool use_blocks = blocklens.size() > 0;
  if (use_blocks) {
    run_id.reserve(T);
    int filled = 0;
    int rid = 0;
    for (int i = 0; i < blocklens.size(); ++i) {
      int len = blocklens[i];
      if (len <= 0) stop("blocklens must contain positive run lengths.");
      for (int t = 0; t < len && filled < T; ++t) {
        run_id.push_back(rid);
        ++filled;
      }
      ++rid;
    }
    if (filled != T) stop("Sum of blocklens must match the number of rows in U.");
  }

  NumericVector v(r, 0.0);
  std::vector< NumericMatrix > accumulators;
  accumulators.reserve(K * L);
  for (int idx = 0; idx < K * L; ++idx) {
    NumericMatrix M(r, r);
    std::fill(M.begin(), M.end(), 0.0);
    accumulators.emplace_back(M);
  }

  auto y_at = [&](int t, int off, int j) -> double {
    int idx = t - off;
    if (idx < 0) idx = 0;
    if (idx >= T) idx = T - 1;
    return U(idx, j);
  };

  std::vector< std::vector<double> > lagged_cache(L, std::vector<double>(r));
  std::vector<double> ut(r);

  for (int t = 0; t < T; ++t) {
    for (int j = 0; j < r; ++j) {
      ut[j] = U(t, j);
    }

    for (int i = 0; i < r; ++i) {
      double ui = ut[i];
      v[i] = a * v[i] + one_minus_a * (ui * ui);
    }

    for (int li = 0; li < L; ++li) {
      const int off = offsets[li];
      for (int j = 0; j < r; ++j) {
        lagged_cache[li][j] = y_at(t, off, j);
      }
    }

    for (int k = 0; k < K; ++k) {
      double g = PK(t, k);
      for (int li = 0; li < L; ++li) {
        NumericMatrix &A = accumulators[k * L + li];
        for (int i = 0; i < r; ++i) {
          for (int j = 0; j < r; ++j) {
            A(i, j) *= a;
          }
        }
        if (R_IsNA(g)) continue;

        int idx = t - offsets[li];
        if (idx < 0) idx = 0;
        if (idx >= T) idx = T - 1;
        if (use_blocks && run_id[t] != run_id[idx]) continue;

        const std::vector<double> &y = lagged_cache[li];
        const double weight = one_minus_a * g;
        for (int i = 0; i < r; ++i) {
          double wi = weight * ut[i];
          for (int j = 0; j < r; ++j) {
            A(i, j) += wi * y[j];
          }
        }
      }
    }
  }

  if (normalized) {
    for (int idx = 0; idx < K * L; ++idx) {
      NumericMatrix &A = accumulators[idx];
      for (int i = 0; i < r; ++i) {
        double vi = std::max(v[i], 1e-12);
        for (int j = 0; j < r; ++j) {
          double denom = std::sqrt(vi * std::max(v[j], 1e-12));
          if (denom > 0) {
            A(i, j) /= denom;
          } else {
            A(i, j) = NA_REAL;
          }
        }
      }
    }
  }

  List out(K * L);
  for (int k = 0; k < K; ++k) {
    for (int li = 0; li < L; ++li) {
      out[k * L + li] = accumulators[k * L + li];
    }
  }
  out.attr("K") = K;
  out.attr("L") = L;
  out.attr("offsets") = offsets;
  return out;
}

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
    if (R_IsNA(xv) || R_IsNA(yv)) {
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
      double denom = std::sqrt(std::max(Sxx * Syy, 1e-12));
      out[t] = denom > 0 ? (Sxy / denom) : NA_REAL;
    }
  }

  return out;
}
