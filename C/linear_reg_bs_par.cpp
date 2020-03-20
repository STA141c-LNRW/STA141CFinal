#include <RcppArmadillo.h>
#include <boost/math/distributions/students_t.hpp>

using namespace Rcpp;

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]

/*
library(Rcpp)
library(RcppArmadillo)
x = data.frame(1:3000 + runif(3000,0,100), 3001:6000 + runif(3000,0,100))
y = 6001:9000 + runif(3000, 0, 100)
s = 10
r = 100
 sourceCpp("C/linear_reg_bs_par.cpp")
 linear_reg_bs_par_C(x,y,s,r)
*/

// [[Rcpp::export]]
void linear_reg_bs_par_C(DataFrame x, arma::colvec y, int s, int r){
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
  for (int i = 0; i < s; i++){
    for (int j = 0; j < keepsplit; j++){
      for(int k = 0; k < p; k++){
        // x_samples(split, p ,s)
        // y_samples(split, s)
        Rcout << x_samples(j,k,i) << ",";
      }
      Rcout << ",y = " << y_samples(j,i);
      Rcout << std::endl;
    }
  }
  // x_samples(split, p, s)
  // x_samples(split, p, i)
  // need to count split
  for (int i = 0; i < s; i++){
    std::vector<int> f(r);
    int n_sub = keepsplit;
    if (x_samples(keepsplit, 1, i) != 1){
      n_sub = n_sub - 1;
    }
    std::fill(f.begin(), f.end(), 1);
    std::discrete_distribution<int> distribution(f.begin(), f.end());
    std::default_random_engine generator(i);
    std::fill(f.begin(), f.end(), 0);
    for (int j = 0; j < n; j++){
      int number = distribution(generator);
      ++f[number];
    }
    for (freq:f){

    }
  }
}