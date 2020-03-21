#include <RcppArmadillo.h>
using namespace Rcpp;

//' This function takes in a list of linear regression coefficient estimates generated
//' by a Bag of Little Bootstraps procedure, and a dataframe of observations without the
//' response variable. The response variable for each observation is predicted using
//' each vector of coefficient estimates for each sample. Then, empirical prediction
//' intervals and point estimates for the response variable of each observation are
//' determined for each sample. Afterwards, the endpoints of all intervals are averaged
//' to form overall prediction intervals, and point estimates are averaged to form
//' overall predictions. It should be noted that the prediction intervals are not
//' multiple prediction intervals. For Bonferroni-corrected prediction intervals, divide
//' the desired value of alpha by the number of observations. The difference between
//' this function and PI is that this function is written in C++ instead of R for faster
//' performance.
//'
//' @param lrbs A linear_reg_bs or linear_reg_bs_par object containing BLB regression
//' coefficient estimates.
//' @param x A dataframe of the explanatory variables of unseen observations.
//' @param alpha The significance level. Default value is 0.05.
//' @return The prediction intervals and estimates for the response variable of each
//' unseen observation.
//' @export
//'
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix PI_C(List lrbs, DataFrame x, double alpha){
  double lowerq = alpha/2;
  double upperq = 1 - alpha/2;
  List coefs = lrbs["bootstrap_coefficient_estimates"];
  NumericMatrix temp = coefs[0];
  int replications = temp.ncol();
  int n = x.nrows();
  int p = x.size() + 1;
  arma::mat z(n,p);
  std::fill(z.begin(), z.begin() + n, 1);
  for (int i = 1; i < p; i++){
    NumericVector temp = x[i-1];
    std::copy(temp.begin(), temp.end(), z.begin() + i*n);
  }
  arma::cube preds(replications, coefs.size(), x.nrows(), arma::fill::zeros);
  NumericMatrix prediction_intervals(x.nrows(), 3);
  for (int i = 0; i < coefs.size(); i++){
    NumericMatrix sub_coefs = coefs[i];
    for (int j = 0; j < replications; j++){
      for (int k = 0; k < n; k++){
        for (int m = 0; m < p; m++){
          preds(j,i,k) += z(k,m)*sub_coefs(m,j);
        }
      }
    }
  }
  for (int i = 0; i < x.nrows(); i++){
    double lowerE = 0.0;
    double estimatesE = 0.0;
    double upperE = 0.0;
    for (int k = 0; k < coefs.size(); k++){
        double sum = 0;
        std::vector<double> subset(replications);
        for (int j = 0; j < replications; j++){
          subset[j] = preds(j,k,i);
          sum += subset[j];
        }
        std::sort(subset.begin(), subset.end());
        lowerE += subset[subset.size()*lowerq];
        upperE +=  subset[subset.size()*upperq];
        estimatesE += sum/subset.size();
    }
    prediction_intervals(i, 0) = lowerE/coefs.size();
    prediction_intervals(i, 1) = estimatesE/coefs.size();
    prediction_intervals(i, 2) = upperE/coefs.size();
  }
  colnames(prediction_intervals) = CharacterVector(
    {"Lower_Bounds", "Estimates", "Upper_Bounds"});
  return prediction_intervals;
}

