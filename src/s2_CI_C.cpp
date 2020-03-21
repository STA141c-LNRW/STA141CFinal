#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
NumericVector s2_CI_C(List lrbs, double alpha){
  double lowerq = alpha/2;
  double upperq = 1 - alpha/2;
  List coefs = lrbs["bootstrap_s2_estimates"];
  double lower = 0.0;
  double upper = 0.0;
  double estimates = 0.0;
  for (int i = 0; i < coefs.size(); i++){
    NumericVector subset = coefs[i];
    double sum = 0;
    for (int i = 0; i < subset.size(); i++){
      sum += subset[i];
    }
    std::sort(subset.begin(), subset.end());
    lower += subset[subset.size()*lowerq];
    upper +=  subset[subset.size()*upperq];
    estimates += sum/subset.size();
  }
  lower /= coefs.size();
  estimates /= coefs.size();
  upper /= coefs.size();
  NumericVector prediction_intervals = NumericVector::create(
    _["Lower_Bound"] = lower,
    _["Estimate"] = estimates,
    _["Upper_Bound"] = upper
  );
  return prediction_intervals;
}
