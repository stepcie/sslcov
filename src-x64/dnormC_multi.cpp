#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
const double log2pi2 = log(2.0 * M_PI)/2;

//' C++ implementation of univariate Normal of mean zero pdf for multiple inputs
//'
//'@param x data matrix of dimension p x n, p being the dimension of the
//'data and n the number of data points
//'@param sd the standard deviation
//'@param logical flag for returning the log of the probability density
//'function. Defaults is \code{TRUE}
//'@return matrix of densities of dimension \code{p x n}.
//'@export
// [[Rcpp::export]]
NumericMatrix dnormC_multi(NumericVector x,
                           NumericVector m,
                           double sd,
                           bool Log=false){

  vec xx = as<vec>(x);
  vec mm = as<vec>(m);
  int p = xx.n_elem;
  int n = mm.n_elem;
  NumericMatrix y(p,n);
  double constant = - log(sd) - log2pi2;

  for (int i=0; i < p; i++) {
    for (int j=0; j < n; j++) {
      if (!Log) {
        y(i,j) = exp(-0.5*pow((xx(i)-mm(j))/sd,2) + constant);
      } else{
        y(i,j) = -0.5*pow((xx(i)-mm(j))/sd,2)  + constant;
      }
    }
  }

  return y;

}
