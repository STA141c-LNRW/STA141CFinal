
#include <RcppArmadillo.h>
#include <boost/math/distributions/students_t.hpp>

using namespace Rcpp;

// [[Rcpp::depends(BH)]]
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double tc(double quant, int tf){
  boost::math::students_t dist(tf);
  return quantile(dist, quant);
}

// [[Rcpp::export]]
SEXP tr(double quant, int tf){
  Function qt("qt");
  return qt(quant, tf);
}

// [[Rcpp::export]]
NumericVector lr_coefficient_CI_C(List lr, double alpha){
  int n = lr["n"];
  int p = lr["p"];
  double q = tc(1 - alpha/2, n - p);
  NumericMatrix inverse = lr["xtx_inverse"];
  double sigma = lr["sigma_2_hat"];
  NumericVector coefficients = lr["coefficients"];
  NumericMatrix final(coefficients.size(), 2);
  for (int i = 0; i < coefficients.size(); i++){
    double var_coeffs = sigma * inverse(i, i);
    final[i] = coefficients[i] - q*sqrt(var_coeffs);
    final[i+coefficients.size()] = coefficients[i] + q*sqrt(var_coeffs);
  }
  colnames(final) = CharacterVector::create("Lower","Upper");
  return final;
}

// [[Rcpp::export]]
double calc_slopeC(DataFrame df){
  NumericVector x = df["x"];
  int length = x.size();
  NumericVector y = df["y"];
  double meanX = 0;
  double meanY = 0;
  for (int i = 0; i < length; i++){
    meanX += x[i];
    meanY += y[i];
  }
  meanX /= length;
  meanY /= length;
  double sumSXX = 0;
  double sumSXY = 0;
  for (int i = 0; i < length; i++){
    double temp = x[i] - meanX;
    sumSXX += temp*temp;
    sumSXY += temp*(y[i] - meanY);
  }
  return sumSXY/sumSXX;
}

// [[Rcpp::export]]
NumericVector lr_prediction_interval_C(List lr, arma::colvec x, double alpha){
  x.insert_rows(0, 1);
  x[0] = 1;
  double sigma = lr["sigma_2_hat"];
  NumericMatrix tempinv = lr["xtx_inverse"];
  arma::mat xtx_inverse = arma::mat(tempinv.begin(), tempinv.nrow(),
                                    tempinv.ncol(),
                                    false);
  double var = (1 + trans(x) * xtx_inverse * x)[0] * sigma;
  int n = lr["n"];
  int p = lr["p"];
  double qt = tc(1 - alpha/2, n - p);
  NumericMatrix temp = lr["coefficients"];
  double mean_res = 0;
  for (unsigned int i = 0; i < x.n_rows; i++){
    mean_res += temp[i]*x[i];
  }
  double lower = mean_res - (qt * sqrt(var));
  double upper = mean_res + (qt * sqrt(var));
  NumericVector final = NumericVector::create(lower,mean_res,upper);
  final.attr("dim") = Dimension(1,3);
  return final;
}

// [[Rcpp::export]]
NumericVector lr_prediction_interval_C2(List lr, NumericVector x, double alpha){
  x.insert(0, 1);
  double sigma = lr["sigma_2_hat"];
  NumericMatrix tempinv = lr["xtx_inverse"];
  NumericVector y(tempinv.ncol());
  for (int i = 0; i < tempinv.ncol(); i++){
    y[i] = std::inner_product(x.begin(), x.end(),
                       tempinv.begin()+i*tempinv.ncol(), 0.);
  }
  double t = std::inner_product(y.begin(), y.end(),
                                   x.begin(), 0.);
  double var = (1 + t) * sigma;
  int n = lr["n"];
  int p = lr["p"];
  double qt = tc(1 - alpha/2, n - p);
  NumericMatrix temp = lr["coefficients"];
  double mean_res = 0;
  for (unsigned int i = 0; i < x.size(); i++){
    mean_res += temp[i]*x[i];
  }
  double lower = mean_res - (qt * sqrt(var));
  double upper = mean_res + (qt * sqrt(var));
  NumericVector final = NumericVector::create(lower,mean_res,upper);
  final.attr("dim") = Dimension(1,3);
  return final;
}


// [[Rcpp::export]]
arma::mat armamultiply(arma::mat x, arma::mat y){
  return x*y;
}


// [[Rcpp::export]]
arma::mat armainverse(arma::mat x){
  return arma::inv(x);
}

// [[Rcpp::export]]
NumericMatrix multiply(NumericMatrix x, NumericMatrix y){
  int one = x.nrow();
  int two = y.ncol();
  int three = x.ncol();
  NumericMatrix z(one,two);
  for (int i = 0; i < one; i++){
    for (int j = 0; j < two; j++){
      for (int k = 0; k < three; k++){
        z(i,j) += x(i,k) * y(k,j);
      }
    }
  }
  return z;
}

// [[Rcpp::export]]
NumericMatrix multiply2(NumericMatrix x, NumericMatrix y){
  int one = x.nrow();
  int two = y.ncol();
  int travel = x.ncol();
  NumericMatrix z(one,two);
  NumericMatrix transx(travel, one);
  for (int i = 0; i < travel; i++){
    for (int j = 0; j < one; j++){
      transx(i, j) = x(j, i);
    }
  }
  for (int i = 0; i < one; i++){
    for (int j = 0; j < two; j++){
      z(i,j) = std::inner_product(transx.begin()+i*travel, transx.begin()+(i+1)*travel, y.begin()+j*travel, 0.);
    }
  }
  return z;
}

// [[Rcpp::export]]
NumericMatrix inverse(NumericMatrix x){
  int i, j, k;
  // Gauss-Jordan Method
  int order = x.nrow();
  NumericMatrix y(order,order*2);
  double temp;
  for (i = 0; i < order; i++){
    for (j = 0; j < order; j++){
      y(i,j) = x(i,j);
    }
  }
  for (i = 0; i < order; i++){
    for (j = order; j < 2*order; j++){
      if (j == (i + order)){
        y(i,j) = 1;
      }
    }
  }
  for (i = order - 1; i > 0; i--)
  {
    if (y(i-1, 1) < y(i, 1)){
      for (j = 0; j < order * 2; j++)
      {
        temp = y(i, j);
        y(i, j) = y(i - 1, j);
        y(i - 1, j) = temp;
      }
    }
  }
  for (i = 0; i < order; i++){
    for (j = 0; j < order; j++){
      if (j != i){
        temp = y(j, i)/y(i, i);
        for (k = 0; k < order*2; k++){
          y(j, k) = y(j, k) - y(i, k)*temp;
        }
      }
    }
  }
  for (i = 0; i < order; i++){
    temp = y(i, i);
    for (j = 0; j < order*2; j++){
      y(i, j) = y(i, j)/temp;
    }
  }
  return y(Range(0,order-1),Range(order,2*order-1));
}


NumericMatrix multiply(NumericMatrix x, NumericMatrix y);
NumericMatrix inverse(NumericMatrix x);
// [[Rcpp::export]]
List linear_regC(NumericMatrix x, NumericVector y){
  int n = x.nrow();
  int p = x.ncol() + 1;
  NumericMatrix m1(n,p);
  std::fill(m1.begin(), m1.begin() + n, 1);
  std::copy(x.begin(), x.end(), m1.begin() + n);
  NumericMatrix tm1 = transpose(m1);
  y.attr("dim") = Dimension(n,1);
  // Solution to OLS is B = (X^T * X)^-1 * X^T * y
  NumericMatrix xtxinverse = inverse(multiply(tm1, m1));
  NumericMatrix coeffs = multiply(multiply(xtxinverse, tm1),
                                  as<NumericMatrix>(y));
  NumericMatrix fv = multiply(m1, coeffs);
  NumericMatrix res(n,1);
  for (int i = 0; i < n; i++){
    res[i] = y[i] - fv[i];
  }
  double s2 = 0;
  for (int i = 0; i < res.size(); i++){
    s2 += res[i]*res[i];
  }
  s2 = s2/(n-p);
  return List::create(_["coefficients"] = coeffs,
                      _["fitted_values"] = fv,
                      _["residuals"] = res,
                      _["sigma_2_hat"] = s2,
                      _["xtx_inverse"] = xtxinverse,
                      _["n"] = n,
                      _["p"] = p);
}

NumericMatrix multiply2(NumericMatrix x, NumericMatrix y);
NumericMatrix inverse(NumericMatrix x);
// [[Rcpp::export]]
List linear_regC2(NumericMatrix x, NumericVector y){
  int n = x.nrow();
  int p = x.ncol() + 1;
  NumericMatrix m1(n,p);
  std::fill(m1.begin(), m1.begin() + n, 1);
  std::copy(x.begin(), x.end(), m1.begin() + n);
  y.attr("dim") = Dimension(n,1);
  NumericMatrix z(p,p);
  NumericMatrix tm1(p, n);
  for (int i = 0; i < tm1.nrow(); i++){
    for (int j = 0; j < tm1.ncol(); j++){
      tm1(i, j) = m1(j, i);
    }
  }
  for (int i = 0; i < p; i++){
    for (int j = 0; j < p; j++){
      z(i,j) = std::inner_product(m1.begin()+i*n, m1.begin()+(i+1)*n, m1.begin()+j*n, 0.);
    }
  }
  // Solution to OLS is B = (X^T * X)^-1 * X^T * y
  // xtxinverse is p * p
  // tm1 is p * n
  // intermed is n * p and we have intermed times y (y is n*1)
  NumericMatrix xtxinverse = inverse(z);
  NumericMatrix intermed = transpose(multiply2(xtxinverse, tm1));
  NumericMatrix coeffs(p, 1);
  for (int i = 0; i < p; i++){
    coeffs[i] = std::inner_product(intermed.begin()+i*intermed.nrow(),
                               intermed.begin()+(i+1)*intermed.nrow(), y.begin(), 0.);
  }
  NumericMatrix fv(n, 1);
  NumericMatrix res(n,1);
  double s2 = 0;
  for (int i = 0; i < n; i++){
    fv[i] = std::inner_product(tm1.begin()+i*p, tm1.begin()+(i+1)*p, coeffs.begin(), 0.);
    res[i] = y[i] - fv[i];
    s2 += res[i]*res[i];
  }
  s2 = s2/(n-p);
  return List::create(_["coefficients"] = coeffs,
                      _["fitted_values"] = fv,
                      _["residuals"] = res,
                      _["sigma_2_hat"] = s2,
                      _["xtx_inverse"] = xtxinverse,
                      _["n"] = n,
                      _["p"] = p);
}

// [[Rcpp::export]]
List linear_regC3(arma::mat x, arma::vec y){
  int n = x.n_rows;
  int p = x.n_cols + 1;
  arma::mat m1(n,p);
  std::fill(m1.begin(), m1.begin() + n, 1);
  std::copy(x.begin(), x.end(), m1.begin() + n);
  arma::mat tm1 = trans(m1);
  arma::mat xtxinverse = arma::inv(tm1* m1);
  arma::mat coeffs = xtxinverse * tm1 * arma::mat(y);
  arma::mat fv = m1 * coeffs;
  NumericMatrix res(n,1);
  for (int i = 0; i < n; i++){
    res[i] = y[i] - fv[i];
  }
  double s2 = 0;
  for (int i = 0; i < res.size(); i++){
    s2 += res[i]*res[i];
  }
  s2 = s2/(n-p);
  return List::create(_["coefficients"] = as<NumericMatrix>(wrap(coeffs)),
                      _["fitted_values"] = as<NumericMatrix>(wrap(fv)),
                      _["residuals"] = res,
                      _["sigma_2_hat"] = s2,
                      _["xtx_inverse"] = as<NumericMatrix>(wrap(xtxinverse)),
                      _["n"] = n,
                      _["p"] = p);
}

