// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>

using arma::mat;
using arma::sp_mat;
using arma::vec;

static sp_mat convmtx_sparse(const vec& h, const int T) {
  const int L = h.n_elem;
  std::vector<arma::uword> ii;
  std::vector<arma::uword> jj;
  std::vector<double> vv;
  ii.reserve(static_cast<std::size_t>(T) * L);
  jj.reserve(static_cast<std::size_t>(T) * L);
  vv.reserve(static_cast<std::size_t>(T) * L);

  for (int t = 0; t < T; ++t) {
    for (int k = 0; k < L; ++k) {
      int col = t - k;
      if (col >= 0) {
        ii.push_back(static_cast<arma::uword>(t));
        jj.push_back(static_cast<arma::uword>(col));
        vv.push_back(h(k));
      }
    }
  }

  const std::size_t nnz = vv.size();
  arma::umat locations(2, nnz);
  arma::vec values(nnz);
  for (std::size_t idx = 0; idx < nnz; ++idx) {
    locations(0, idx) = ii[idx];
    locations(1, idx) = jj[idx];
    values(idx) = vv[idx];
  }

  return sp_mat(locations, values, T, T);
}

static sp_mat diff1_R(const int T) {
  std::vector<arma::uword> ii;
  std::vector<arma::uword> jj;
  std::vector<double> vv;
  ii.reserve(static_cast<std::size_t>(3) * T);
  jj.reserve(static_cast<std::size_t>(3) * T);
  vv.reserve(static_cast<std::size_t>(3) * T);

  for (int i = 0; i < T; ++i) {
    if (i == 0) {
      ii.push_back(i); jj.push_back(i);   vv.push_back(1.0);
      ii.push_back(i); jj.push_back(i+1); vv.push_back(-1.0);
    } else if (i == T - 1) {
      ii.push_back(i); jj.push_back(i-1); vv.push_back(-1.0);
      ii.push_back(i); jj.push_back(i);   vv.push_back(1.0);
    } else {
      ii.push_back(i);   jj.push_back(i-1); vv.push_back(-1.0);
      ii.push_back(i);   jj.push_back(i);   vv.push_back(2.0);
      ii.push_back(i);   jj.push_back(i+1); vv.push_back(-1.0);
    }
  }

  const std::size_t nnz = vv.size();
  arma::umat locations(2, nnz);
  arma::vec values(nnz);
  for (std::size_t idx = 0; idx < nnz; ++idx) {
    locations(0, idx) = ii[idx];
    locations(1, idx) = jj[idx];
    values(idx) = vv[idx];
  }

  return sp_mat(locations, values, T, T);
}

// [[Rcpp::export]]
arma::mat mppi_deconv_map(const arma::mat& Y, const arma::vec& h, const double lambda) {
  const int T = Y.n_rows;
  const int V = Y.n_cols;
  sp_mat H = convmtx_sparse(h, T);
  sp_mat R = diff1_R(T);
  sp_mat A = H.t() * H + lambda * R;
  arma::mat X(T, V, arma::fill::zeros);
  arma::mat A_dense(A);
  arma::mat Rchol = arma::chol(A_dense);
  for (int v = 0; v < V; ++v) {
    arma::vec b = H.t() * Y.col(v);
    arma::vec y = arma::solve(arma::trimatl(Rchol.t()), b);
    arma::vec x = arma::solve(arma::trimatu(Rchol), y);
    X.col(v) = x;
  }
  return X;
}
