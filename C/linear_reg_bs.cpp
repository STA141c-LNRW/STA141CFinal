#include <RcppArmadillo.h>
#include <boost/math/distributions/students_t.hpp>

using namespace Rcpp;

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]

/*
library(Rcpp)
library(RcppArmadillo)
x = data.frame(1:3000, 3001:6000 + runif(3000, 0, 1000))
y = 6001:9000 + runif(3000, 0, 2000)
s = 10
r = 30
 sourceCpp("C/linear_reg_bs_par.cpp")
 linear_reg_bs_par_C(x,y,s,r)
*/

// [[Rcpp::export]]
List linear_reg_bs_C(DataFrame x, arma::colvec y, int s, int r){
  int n = x.nrows();
  int p = x.size() + 1;
  arma::mat z(n,p);
  std::fill(z.begin(), z.begin() + n, 1);
  for (int i = 1; i < p; i++){
    NumericVector temp = x[i-1];
    std::copy(temp.begin(), temp.end(), z.begin() + i*n);
  }
  std::vector<int> sample_indices(n);
  std::vector<int> samples(s);
  std::iota(std::begin(sample_indices), std::end(sample_indices), 0);
  std::iota(std::begin(samples), std::end(samples), 0);
  std::random_shuffle(sample_indices.begin(), sample_indices.end());
  std::random_shuffle(samples.begin(), samples.end());
  int split = n/s + (n % s != 0);
  int keepsplit = split;
  arma::cube x_samples(split,p,s);
  arma::mat y_samples(split,s);
  bool lower = false;
  if (fmod(n + 0.0, s + 0.0) < 0.0001){
    lower = true;
  }
  int remaining = n;
  int sum = 0;
  for (int i = 0; i < s; i++){
    if (!lower && fmod(remaining + 0.0, s - i + 0.0) < 0.0001){
      lower = true;
      split = split - 1;
    }
    remaining = remaining - split;
    int sindex = samples[i];
    for (int j = 0; j < split; j++){
      for (int k = 0; k < p; k++){
        x_samples(j, k, sindex) = z(sample_indices[sum], k);
      }
      y_samples(j, sindex) = y[sample_indices[sum]];
      sum++;
    }
  }
  // x_samples(split, p, s)
  // x_samples(split, p, i)
  // need to count split
  List bs_coefs(s);
  List bs_s2(s);
  for (int i = 0; i < s; i++){
    NumericMatrix sample_coefs(p,r);
    NumericVector sample_s2(r);
    int n_sub = keepsplit;
    if (x_samples(keepsplit - 1, 1, i) != 1){
       n_sub = n_sub - 1;
    }
    for (int j = 0; j < r; j++){
      std::vector<int> freqs(n_sub);
      std::fill(freqs.begin(), freqs.end(), 1);
      std::discrete_distribution<int> distribution(freqs.begin(), freqs.end());
      std::default_random_engine generator(j);
      std::fill(freqs.begin(), freqs.end(), 0);
      for (int k = 0; k < n; k++){
        int number = distribution(generator);
        ++freqs[number];
      }
      arma::mat x_resamp(n,p);
      arma::mat y_resamp(n,1);
      int row = 0;
      for (int k = 0; k < n_sub; k++){
        int repeat = freqs[k];
        for (int l = 0; l < repeat; l++){
          for (int m = 0; m < p; m++){
            x_resamp(row, m) = x_samples(k, m, i);
          }
          y_resamp(row, 0) = y_samples(k, i);
          row++;
        }
      }
      arma::mat transx_resamp = trans(x_resamp);
      arma::mat xtxinverse = arma::inv(transx_resamp * x_resamp);
      arma::mat coefs = xtxinverse * transx_resamp * y_resamp;
      arma::mat fv = x_resamp * coefs;
      arma::mat res = y_resamp - fv;
      sample_s2(j) = std::inner_product(res.begin(), res.end(), res.begin(), 0.)/
        (n_sub - p);
      for (int k = 0; k < p; k++){
        sample_coefs(k,j) = coefs[k];
      }
    }
    bs_coefs(i) = sample_coefs;
    bs_s2(i) = sample_s2;
  }
  return List::create(_["bootstrap_coefficient_estimates"] = bs_coefs,
                      _["bootstrap_s2_estimates"] = bs_s2);
}
