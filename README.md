# STA141CFinal
Implementation of Multiple Linear Regression using BLB

Final Project for STA 141C at UC Davis, Winter 2020

William Shih, Ricardo Simpao, Nilay Varshney, Luke Yee

This R package provides functions for fitting linear regression models on datasets with continuous response variables through the use of Bag of Little Bootstraps. Users can determine regression coefficients, estimate the variance of errors, and predict new data. Users also have the option to use parallel computing.

The dataset used in the user guide can be found at the link below. Note that the three rightmost columns are to be disregarded as the dataset has four response variables.

archive.ics.uci.edu/ml/datasets/SGEMM+GPU+kernel+performance

# TestRandC.pdf 
file shows the results of benchmarking basic R and C functions to see what to use for the final package.

# User_Guide.pdf 
file shows the user guide on how to use the package.

# package_demo.pdf 
benchmarks the 8 R functions in the package versus the 4 C++ functions in the package.

There are 8 R files corresponding to 8 R functions:

PI.R

coef_CI.R

linear_reg_bs.R

s2_CI.R

PI_par.R

coef_CI_par.R

linear_reg_bs_par.R

s2_CI_Par.R

There are 4 C++ files corresponding to 4 C++ functions:

PI_C.cpp

coef_CI_C.cpp

linear_reg_bs_C.cpp

s2_CI_C.cpp

