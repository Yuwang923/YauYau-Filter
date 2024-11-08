#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
List KolmogorovEW(int Dim, double Dt, arma::vec s, 
                  Function f, Function df, Function h) {
  double Ds = s(1) - s(0);
  int Ns = s.n_elem;
  vec e = ones<vec>(Ns);
  sp_mat K = 0.5 / Ds * spdiags(join_horiz(-e, e), {-1, 1}, Ns, Ns);
  
  // Create the grid
  field<mat> grids(Dim);
  for (int i = 0; i < Dim; ++i) {
    grids(i) = repmat(s, 1, pow(Ns, Dim - i - 1));
  }
  mat s_comb = vectorise(join_horiz(grids));
  
  // Create derivative operators for each dimension
  field<sp_mat> D(Dim);
  for (int i = 0; i < Dim; ++i) {
    D(i) = kron(sp_mat(eye<mat>(static_cast<uword>(pow(Ns, Dim - i - 1)), static_cast<uword>(pow(Ns, Dim - i - 1)))), kron(K, sp_mat(eye<mat>(static_cast<uword>(pow(Ns, i)), static_cast<uword>(pow(Ns, i))))));
  }
  
  mat f_val = as<mat>(f(s_comb));
  vec tmp = -vectorise(f_val);
  mat df_val = as<mat>(df(s_comb));
  mat h_val = as<mat>(h(s_comb));
  vec diag_elements = -(sum(df_val, 1) + 0.5 * sum(square(h_val), 1));
  sp_mat Q = spdiags(diag_elements, 0, pow(Ns, Dim), pow(Ns, Dim));
  for (int i = 0; i < Dim; ++i) {
    Q += spdiags(tmp, 0, pow(Ns, Dim), pow(Ns, Dim)) * D(i);
  }
  sp_mat B = sp_mat(eye<mat>(static_cast<uword>(pow(Ns, Dim)), static_cast<uword>(pow(Ns, Dim)))) + Dt * Q;
  
  // Calculate Lambda
  vec d = -4 * square(sin(linspace<vec>(1, Ns, Ns) * M_PI / (2 * Ns + 2)));
  sp_mat D_kron(static_cast<uword>(pow(Ns, Dim)), static_cast<uword>(pow(Ns, Dim)));
  D_kron.zeros();
  for (int i = 0; i < Dim; ++i) {
    D_kron += kron(sp_mat(eye<mat>(static_cast<uword>(pow(Ns, Dim - i - 1)), static_cast<uword>(pow(Ns, Dim - i - 1)))), kron(spdiags(d, 0, Ns, Ns), sp_mat(eye<mat>(static_cast<uword>(pow(Ns, i)), static_cast<uword>(pow(Ns, i))))));
  }
  vec Lambda = 1 - 0.5 * Dt / (Ds * Ds) * D_kron.diag();
  
  return List::create(Named("Lambda") = Lambda,
                      Named("B") = B,
                      Named("s") = s_comb,
                      Named("Ns") = Ns);
}