#include <Rcpp.h>
using namespace Rcpp;
//'C++ implementation of empirical cdf
//'
//'@param sobs the observation for which the cdf is to be evaluated
//'
//'@param sref the reference observations from which the cdf is to be constructed
//'
//'@useDynLib sslcov
//'@examples
//'library(microbenchmark)
//'fi <- rnorm(n=1000)
//'n_new <- 100
//'fnew <- rnorm(n=n_new) 
//'microbenchmark(#sum.I(c(fi, fnew), FUN=">=", fnew)/n_new,
//'    ecdf(fnew)(c(fi, fnew)),
//'    ecdf_cpp(c(fi, fnew), sort(fnew))/n_new
//')
//'@export
//[[Rcpp::export]]
NumericVector ecdf_cpp(NumericVector sobs, NumericVector sref) {
  int nobs = sobs.size();
  NumericVector ans(nobs);
  std::sort(sref.begin(), sref.end());
  for (int i = 0; i < nobs; ++i){
    ans[i] = (std::upper_bound(sref.begin(), sref.end(), sobs[i]) - sref.begin());
  }
  return ans/((double) nobs);
}
