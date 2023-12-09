#include <iostream>
#include <vector>
#include <random>
#include <Rcpp.h>
using namespace Rcpp;

//' @title A ARMA series sampler using Rcpp
//' @description A ARMA series sampler using Rcpp
//' @param p  order of autoregression
//' @param q  order of moving average
//' @param n  length of time series
//' @param arParams a vector of autoregression parameters
//' @param maParams a vector of moving average parameters
//' @return the desired arma series
//' @useDynLib SA23204170
//' @examples
//' \dontrun{
//' error<- generateARMA(3,4,1000,c(.5,-.2),.3)
//' }
//' @export
// [[Rcpp::export]]
NumericVector generateARMA(int p, int q, int n, const Rcpp::NumericVector& arParams, const Rcpp::NumericVector& maParams) {
  std::vector<double> series(n, 0.0);
  // Seed for random number generation
  std::random_device rd;
  std::mt19937 gen(rd());
  std::normal_distribution<double> distribution(0.0, 1.0);
  // Generating ARMA time series
  for (int i = std::max(p, q); i < n; ++i) {
    double arComp = 0.0, maComp = 0.0;
    // Autoregressive component
    for (int j = 0; j < p; ++j) {
      arComp += arParams[j] * series[i - j - 1];
    }
    // Moving average component
    for (int k = 0; k < q; ++k) {
      maComp += maParams[k] * distribution(gen);
    }
    series[i] = arComp + maComp;
  }
  return Rcpp::wrap(series);
}

//' @title A Gibbs sampler
//' @description A Gibbs sampler
//' @param N number of samples
//' @param n number of trials 
//' @param a non-negative parameters of the Beta distribution.
//' @param b non-negative parameters of the Beta distribution.
//' @useDynLib SA23204170
//' @export
// [[Rcpp::export]]
NumericMatrix gibbsC(int N, int n, int a, int b) {
  NumericMatrix mat(N, 2);
  double x = 0, y = 0;
  for(int i = 0; i < N; i++) {
    x = rbinom(1, 10, y)[0];
    y = rbeta(1, x+a, n-x+b)[0];
    mat(i, 0) = x;
    mat(i, 1) = y;
  }
  return(mat);
}
