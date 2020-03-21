#include <RcppArmadillo.h>
using namespace Rcpp;

//' This function takes in a list of linear regression coefficient estimates generated
//' by a Bag of Little Bootstraps procedure. Then, empirical confidence intervals and
//' point estimates of each coefficient are determined for each subsample. Afterwards,
//' the endpoints of all confidence intervals are averaged to form overall confidence
//' intervals, and point estimates are averaged to form overall estimates. It should be
//' noted that the confidence intervals are not multiple confidence intervals. For
//' Bonferroni-corrected confidence intervals, divide the desired value of alpha by the
//' number of regression coefficients. The difference between this function and coef_CI
//' is that this function is written in C++ instead of R for faster performance.
//'
//' @param lrbs A linear_reg_bs or linear_reg_bs_par object containing BLB regression
//' coefficient estimates.
//' @param alpha The significance level. Default value is 0.05.
//' @return The overall confidence interval for each regression coefficient, along with
//' its overall estimate.
//' @export
//'
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
