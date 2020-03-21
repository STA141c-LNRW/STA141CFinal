#include <RcppArmadillo.h>

using namespace Rcpp;

// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
NumericMatrix coef_CI_C(List lrbs, double alpha){
  double lowerq = alpha/2;
  double upperq = 1 - alpha/2;
  List coefs = lrbs["bootstrap_coefficient_estimates"];
  NumericMatrix one = coefs[0];

  NumericMatrix prediction_intervals(one.nrow(),3);
  std::vector<double> lowerS(one.nrow(), 0);
  std::vector<double> upperS(one.nrow(), 0);
  std::vector<double> sumS(one.nrow(), 0);
  for (int i = 0; i < coefs.size(); i++){
    NumericMatrix s = coefs[i];
    NumericMatrix subset = transpose(s);
    for (int j = 0; j < subset.ncol(); j++){
      std::sort(subset.begin()+j*subset.nrow(),
                subset.begin()+(j+1)*subset.nrow());
      lowerS[j] += subset[j*subset.nrow()+subset.nrow()*lowerq];
      upperS[j] += subset[j*subset.nrow()+subset.nrow()*upperq];
      for (int k = 0; k < subset.nrow(); k++){
        sumS[j] += subset(k, j);
      }
    }
  }
  for (int i = 0; i < one.nrow(); i++){
    prediction_intervals(i, 0) = lowerS[i]/coefs.size();
    prediction_intervals(i, 2) = upperS[i]/coefs.size();
    prediction_intervals(i, 1) = sumS[i]/(coefs.size()*one.ncol());
  }
  rownames(prediction_intervals) = rownames(one);
  colnames(prediction_intervals) = CharacterVector(
  {"Lower_Bounds", "Estimates", "Upper_Bounds"});
  return prediction_intervals;
}
