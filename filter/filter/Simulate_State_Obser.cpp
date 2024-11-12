#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
List Simulate_State_Obser(double Dt, int Ntau, int NtNtau, Function f, Function h, int Dim, Nullable<int> seed = R_NilValue) {
  mat x(NtNtau, Dim, fill::zeros);
  mat y_Dt(NtNtau, Dim, fill::zeros);
  mat y_tau(Ntau, Dim, fill::zeros);
  double sqrtdT = sqrt(Dt);
  
  // Set random seed for reproducibility, use user-provided seed if available
  unsigned int seed_value = seed.isNull() ? std::chrono::system_clock::now().time_since_epoch().count() : as<int>(seed);
  Environment base_env("package:base");
  Function set_seed = base_env["set.seed"];
  set_seed(seed_value);
  Rcout << "Random seed used: " << seed_value << "\n";
  
  // Simulate x
  for (int t = 1; t < NtNtau; ++t) {
    rowvec x_prev = x.row(t - 1);
    NumericVector f_value = as<NumericVector>(f(x_prev));
    rowvec noise = randn<rowvec>(Dim);
    x.row(t) = x_prev + as<rowvec>(f_value) * Dt + sqrtdT * noise;
  }
  
  // Simulate y_Dt
  for (int t = 1; t < NtNtau; ++t) {
    rowvec x_prev = x.row(t);
    NumericVector h_value = as<NumericVector>(h(x_prev));
    rowvec noise = randn<rowvec>(Dim);
    y_Dt.row(t) = y_Dt.row(t - 1) + as<rowvec>(h_value) * Dt + sqrtdT * noise;
  }
  
  // Simulate y_tau
  for (int n = 1; n < Ntau; ++n) {
    int t = n * (NtNtau / Ntau);
    y_tau.row(n) = y_Dt.row(t);
  }
  
  return List::create(Named("x") = x,
                      Named("y_Dt") = y_Dt,
                      Named("y_tau") = y_tau,
                      Named("seed") = seed_value);
}