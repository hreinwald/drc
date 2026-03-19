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

Christian Ritz

## Examples

``` r
## Fitting a four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
vcov(ryegrass.m1)
#>              [,1]        [,2]        [,3]         [,4]
#> [1,]  0.216282967  0.04601511 -0.03504683 -0.003763692
#> [2,]  0.046015113  0.04502563 -0.00471192 -0.016918440
#> [3,] -0.035046835 -0.00471192  0.03555759 -0.012868772
#> [4,] -0.003763692 -0.01691844 -0.01286877  0.034496126
vcov(ryegrass.m1, corr = TRUE)
#>             [,1]       [,2]       [,3]        [,4]
#> [1,]  1.00000000  0.4662936 -0.3996423 -0.04357304
#> [2,]  0.46629357  1.0000000 -0.1177611 -0.42928455
#> [3,] -0.39964231 -0.1177611  1.0000000 -0.36743943
#> [4,] -0.04357304 -0.4292845 -0.3674394  1.00000000
```
