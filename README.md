
# SoftHSMOTRBART

<!-- badges: start -->
<!-- badges: end -->

This package contains an implementation of Bayesian Additive Regression Trees with Model Trees (Prado et al. 2022) and a new adjusted variablion of Hierarchical Shrinkage (Agrawal 2022). This package accounts for the Hierarchical Shrinkage aspect of models within MCMC (drawing parameters, calculating marginal likelihoods etc.), rather than applying shrinkage post-hoc as in Agrawal et al. (2022). The new variation combines soft splitting (Linero and Yang 2018) with hierarchical shrinkage.

[Agarwal, A., Tan, Y. S., Ronen, O., Singh, C., & Yu, B. (2022). _Hierarchical Shrinkage: improving the accuracy and interpretability of tree-based methods_ .](https://arxiv.org/abs/2202.00858)


[Linero, Antonio R., and Yun Yang. _Bayesian regression tree ensembles that adapt to smoothness and sparsity_. Journal of the Royal Statistical Society: Series B (Statistical Methodology) 80.5 (2018): 1087-1110.](https://doi.org/10.1111/rssb.12293)

[Prado, E.B., Moral, R.A. & Parnell, A.C. _Bayesian additive regression trees with model trees_. Statistics and Computing 31, 20 (2021)](https://mural.maynoothuniversity.ie/15498/1/Prado2021_Article_BayesianAdditiveRegressionTree.pdf)


## Installation

You can install the development version of SoftHSMOTRBART like so:

``` r
library(devtools)
install_github("EoghanONeill/SoftHSMOTRBART")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(SoftHSMOTRBART)
## basic example code

library(mvtnorm)
library(stats)
library(MCMCpack)
library(truncnorm)

library(profvis)

# Simulate a Friedman data set
friedman_data = function(n, num_cov, sd_error){
  x = matrix(runif(n*num_cov),n,num_cov)
  y = 10*sin(pi*x[,1]*x[,2]) + 20*(x[,3]-.5)^2+10*x[,4]+5*x[,5] + rnorm(n, sd=sd_error)
  return(list(y = y,
              x = x))
}
# Training data
data = friedman_data(100, 10, 1)
y = data$y
x = data$x

# Test data
data_test = friedman_data(100, 10, 1)
y.test = data_test$y
x.test = data_test$x

# Run MOTR-BART
set.seed(55)
# profvis({
fit.motr.bart.h = hs_motr_bart(x, y, ntrees = 10, nburn = 10, npost = 10)
# })

# save(fit.motr.bart, file = "smotr_result.Rdata")

# test_results = test_function(x.test, fit.motr.bart.h)


# cor(y.test, rowMeans(test_results$predictions))

```

