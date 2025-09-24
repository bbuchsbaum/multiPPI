// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::plugins(openmp)]]
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif
using namespace Rcpp;
using namespace arma;

// Orthonormal DCT-II basis (T x K)
static mat dct_basis(const int T, const int K) {
  mat X(T, K, fill::zeros);
  const double scale0 = std::sqrt(1.0 / T);
  const double scale  = std::sqrt(2.0 / T);
  for (int t = 0; t < T; ++t) {
    for (int k = 0; k < K; ++k) {
      double s = (k == 0) ? scale0 : scale;
      X(t, k) = s * std::cos(M_PI * (t + 0.5) * k / T);
    }
  }
  return X;
}

// Convolve each column of A (length T) with kernel h (length L), 'same' length T
static mat conv_cols_same(const mat &A, const vec &h) {
  const int T = A.n_rows, K = A.n_cols;
  const int L = h.n_elem;
  mat out(T, K, fill::zeros);
  for (int k = 0; k < K; ++k) {
    for (int t = 0; t < T; ++t) {
      double acc = 0.0;
      // y[t] = sum_{tau=0..L-1} h[tau] * a[t - tau], with a[u]=0 for u<0
      int max_tau = std::min(L-1, t);
      for (int tau = 0; tau <= max_tau; ++tau) {
        acc += h(tau) * A(t - tau, k);
      }
      out(t, k) = acc;
    }
  }
  return out;
}

// Compute per-column GCV-optimal lambda on a log-grid using eigen trick
static double gcv_select_lambda(const vec &d, const vec &c, const double yy,
                                const int T, const int ngrid,
                                const double lam_lo, const double lam_hi) {
  double best_lam = lam_lo;
  double best_gcv = std::numeric_limits<double>::infinity();
  for (int i = 0; i < ngrid; ++i) {
    double t = (double)i / (double)(ngrid - 1);
    double lam = std::exp(std::log(lam_lo) * (1 - t) + std::log(lam_hi) * t);
    // rss = ||y||^2 - 2 F^T beta + beta^T G beta
    // with eigen trick: c = U' Q^{-1/2} F, d = eig(A), A = Q^{-1/2} G Q^{-1/2}
    // F^T beta = sum c_i^2 / (d_i + lam)
    // beta^T G beta = sum c_i^2 * d_i / (d_i + lam)^2
    double Ft_b = 0.0, bGb = 0.0, trS = 0.0;
    for (uword j = 0; j < d.n_elem; ++j) {
      double denom = d(j) + lam;
      double ci2 = c(j) * c(j);
      Ft_b += ci2 / denom;
      bGb += ci2 * d(j) / (denom * denom);
      trS  += d(j) / denom;
    }
    double rss = yy - 2.0 * Ft_b + bGb;
    double denom_gcv = (double)T - trS;
    double gcv = rss / (denom_gcv * denom_gcv + 1e-12);
    if (gcv < best_gcv) { best_gcv = gcv; best_lam = lam; }
  }
  return best_lam;
}

// [[Rcpp::export]]
Rcpp::List deconv_dct_multi(const arma::mat &Y,     // T x R signals (prewhitened, or basis-projected)
                            const arma::vec &h,     // HRF kernel (L, in TR units)
                            const int K = 64,       // DCT components
                            Rcpp::Nullable<Rcpp::NumericVector> q_in = R_NilValue,
                            const std::string method = "gcv", // "gcv" or "fixed"
                            const double lambda_fixed = 1e-1,
                            const int ngrid = 30,
                            const double lam_lo = 1e-6,
                            const double lam_hi = 1e+6) {
  const int T = Y.n_rows, R = Y.n_cols;
  if (K <= 0 || K > T) stop("K must be in 1..T");
  // DCT basis and its HRF-convolved version
  mat X  = dct_basis(T, K);        // T x K
  mat HX = conv_cols_same(X, h);   // T x K
  // Gram and eigen trick with frequency prior Q (diagonal)
  mat G = HX.t() * HX;             // K x K
  vec q = arma::ones<vec>(K);
  if (q_in.isNotNull()) {
    Rcpp::NumericVector qR(q_in);
    if ((int)qR.size() == K) {
      for (int i = 0; i < K; ++i) q(i) = std::max(1e-12, (double)qR[i]);
    }
  }
  // A = Q^{-1/2} G Q^{-1/2} = (D * G * D), with D = diag(1/sqrt(q))
  vec invsqrtq = 1.0 / arma::sqrt(q);
  mat D = diagmat(invsqrtq);
  mat A = D * G * D;
  // eigen decomposition A = U diag(d) U'
  vec d; mat U;
  eig_sym(d, U, A);
  // Precompute for speed
  mat UtD = U.t() * D;             // K x K
  // Output containers
  mat Beta(K, R, fill::zeros);
  mat Uest(T, R, fill::zeros);
  vec lam_used(R, fill::zeros);
  // For each column (optionally parallelized across OpenMP threads)
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
  for (int r = 0; r < R; ++r) {
    vec y = Y.col(r);
    vec F = HX.t() * y;                      // K x 1
    double yy = dot(y, y);
    // c = U' Q^{-1/2} F
    vec c = UtD * F;
    double lam = lambda_fixed;
    if (method == "gcv") {
      lam = gcv_select_lambda(d, c, yy, T, ngrid, lam_lo, lam_hi);
    }
    // beta = Q^{-1/2} U diag(1/(d + lam)) c
    vec inv = 1.0 / (d + lam);
    vec tmp = inv % c;                       // elementwise
    vec beta = D * (U * tmp);
    Beta.col(r) = beta;
    Uest.col(r) = X * beta;                  // neural estimate in time domain
    lam_used(r) = lam;
  }
  return Rcpp::List::create(
    Rcpp::Named("U")      = Uest,
    Rcpp::Named("beta")   = Beta,
    Rcpp::Named("lambda") = lam_used
  );
}
