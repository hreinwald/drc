# Calculating variance-covariance matrix for objects of class 'drc'

`vcov` returns the estimated variance-covariance matrix for the
parameters in the non-linear function.

## Usage

``` r
# S3 method for class 'drc'
vcov(object, ..., corr = FALSE, od = FALSE, pool = TRUE, unscaled = FALSE)
```

## Arguments

- object:

  an object of class 'drc'.

- ...:

  additional arguments.

- corr:

  logical. If TRUE a correlation matrix is returned.

- od:

  logical. If TRUE adjustment for over-dispersion is used. This argument
  only makes a difference for binomial data.

- pool:

  logical. If TRUE curves are pooled. Otherwise they are not. This
  argument only works for models with independently fitted curves as
  specified in
  [`drm`](https://hreinwald.github.io/drc/reference/drm.md).

- unscaled:

  logical. If TRUE the unscaled variance-covariance is returned. This
  argument only makes a difference for continuous data.

## Value

A matrix of estimated variances and covariances.

## Author

Christian Ritz, Hannes Reinwald

## Examples

``` r
## Fitting a four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
vcov(ryegrass.m1)
#>              [,1]         [,2]         [,3]         [,4]
#> [1,]  0.216283437  0.046015237 -0.035046934 -0.003763696
#> [2,]  0.046015237  0.045025664 -0.004711938 -0.016918446
#> [3,] -0.035046934 -0.004711938  0.035557624 -0.012868789
#> [4,] -0.003763696 -0.016918446 -0.012868789  0.034496143
vcov(ryegrass.m1, corr = TRUE)
#>             [,1]       [,2]       [,3]        [,4]
#> [1,]  1.00000000  0.4662942 -0.3996428 -0.04357302
#> [2,]  0.46629415  1.0000000 -0.1177614 -0.42928444
#> [3,] -0.39964284 -0.1177614  1.0000000 -0.36743966
#> [4,] -0.04357302 -0.4292844 -0.3674397  1.00000000
```
