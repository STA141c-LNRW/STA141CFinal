#include <RcppArmadillo.h>

using namespace Rcpp;

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

