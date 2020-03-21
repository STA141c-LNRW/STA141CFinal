
//[[Rcpp::depends(BH)]]

#include <Rcpp.h>
#include <stdlib.h>
#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <math.h>
#include <cmath>

using namespace Rcpp;

//[[Rcpp::export]]
NumericVector sammy(NumericVector y, int x){
  NumericVector sample;
  for(int i=0; i<x; i++){
    int index=rand() % y.size();
    sample.insert(i,y[index]);
  }
  return(sample);
}


//[[Rcpp::export]]
NumericVector cnorm(int n, double meen, double sdev){
  boost::mt19937 generator(123);
  NumericVector sample;
  for(int i=0; i<n; i++){
    boost::normal_distribution<double> normdist(meen, sdev);
    double draw = normdist(generator);
    sample.insert(i,draw);
  }
  return(sample);
}



//[[Rcpp::export]]
NumericVector cgamma(int n, double alpha, double beta){
  boost::mt19937 generator(123);
  NumericVector sample;
  for(int i=0; i<n; i++){
    boost::gamma_distribution<double> gdist(alpha, beta);
    int draw = gdist(generator);
    sample.insert(i,draw);
  }
  return(sample);
}







