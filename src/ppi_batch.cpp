// Batch computation of multivariate PPI cross-products
#include <RcppArmadillo.h>
using namespace Rcpp;

static arma::mat make_nan_mat(int r) {
  arma::mat M(r, r);
  M.fill(arma::datum::nan);
  return M;
}

// [[Rcpp::export]]
SEXP ppi_batch_cpp(const arma::mat& U,
                   const arma::mat& G,
                   const IntegerVector& lags,
                   const IntegerVector& blocklens = IntegerVector()) {
  const int T = U.n_rows;
  const int r = U.n_cols;
  if (G.n_rows != T) {
    stop("U and G must have the same number of rows.");
  }
  const int K = G.n_cols;
  const int L = lags.size();

  // Precompute block start indices
  std::vector<int> starts;
  if (blocklens.size() > 0) {
    starts.reserve(blocklens.size());
    int pos = 0;
    for (int i = 0; i < blocklens.size(); ++i) {
      int len = blocklens[i];
      if (len <= 0) stop("blocklens must contain positive lengths.");
      if (pos + len > T) stop("Sum of blocklens exceeds data length.");
      starts.push_back(pos);
      pos += len;
    }
    if (pos != T) stop("Sum of blocklens must equal number of rows in U.");
  }

  List cube_list(L);
  arma::mat denom_mat(K, L, arma::fill::zeros);

  for (int li = 0; li < L; ++li) {
    int lag = lags[li];
    arma::cube Mk(r, r, K, arma::fill::zeros);

    for (int k = 0; k < K; ++k) {
      arma::mat M(r, r, arma::fill::zeros);
      double denom = 0.0;
      bool valid = false;

      if (blocklens.size() == 0) {
        if (lag == 0) {
          arma::vec pk = G.col(k);
          denom = arma::dot(pk, pk);
          if (denom > 0) {
            arma::mat Uw = U.each_col() % pk;
            M = U.t() * Uw;
            valid = true;
          }
        } else if (lag > 0) {
          int span = T - lag;
          if (span > 0) {
            arma::vec pk = G.col(k).head(span);
            denom = arma::dot(pk, pk);
            if (denom > 0) {
              arma::mat R0 = U.rows(0, span - 1);
              arma::mat R1 = U.rows(lag, lag + span - 1);
              arma::mat weighted = R1.each_col() % pk;
              M = R0.t() * weighted;
              valid = true;
            }
          }
        } else { // lag < 0
          int shift = -lag;
          int span = T - shift;
          if (span > 0) {
            arma::vec pk = G.col(k).head(span);
            denom = arma::dot(pk, pk);
            if (denom > 0) {
              arma::mat R0 = U.rows(shift, T - 1);
              arma::mat R1 = U.rows(0, span - 1);
              arma::mat weighted = R1.each_col() % pk;
              M = R0.t() * weighted;
              valid = true;
            }
          }
        }
      } else {
        arma::mat accum(r, r, arma::fill::zeros);
        double denom_total = 0.0;
        for (size_t bi = 0; bi < starts.size(); ++bi) {
          int start = starts[bi];
          int len = blocklens[bi];
          if (lag == 0) {
            arma::vec pk = G.col(k).rows(start, start + len - 1);
            double d = arma::dot(pk, pk);
            if (d > 0) {
              arma::mat Ub = U.rows(start, start + len - 1);
              accum += Ub.t() * (Ub.each_col() % pk);
              denom_total += d;
            }
          } else if (lag > 0) {
            if (len > lag) {
              int span = len - lag;
              arma::vec pk = G.col(k).rows(start, start + span - 1);
              double d = arma::dot(pk, pk);
              if (d > 0) {
                arma::mat R0 = U.rows(start, start + span - 1);
                arma::mat R1 = U.rows(start + lag, start + lag + span - 1);
                accum += R0.t() * (R1.each_col() % pk);
                denom_total += d;
              }
            }
          } else { // lag < 0
            int shift = -lag;
            if (len > shift) {
              int span = len - shift;
              arma::vec pk = G.col(k).rows(start, start + span - 1);
              double d = arma::dot(pk, pk);
              if (d > 0) {
                arma::mat R0 = U.rows(start + shift, start + shift + span - 1);
                arma::mat R1 = U.rows(start, start + span - 1);
                accum += R0.t() * (R1.each_col() % pk);
                denom_total += d;
              }
            }
          }
        }
        denom = denom_total;
        if (denom_total > 0) {
          M = accum;
          valid = true;
        }
      }

      if (valid) {
        M /= denom;
      } else {
        M = make_nan_mat(r);
        denom = 0.0;
      }
      Mk.slice(k) = M;
      denom_mat(k, li) = denom;
    }
    cube_list[li] = Mk;
  }

  return List::create(Named("matrices") = cube_list,
                      Named("denominator") = denom_mat);
}
