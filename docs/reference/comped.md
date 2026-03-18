# Comparison of effective dose values

Comparison of a pair of effective dose values from independent
experiments where only the estimates and their standard errors are
reported.

## Usage

``` r
comped(
  est,
  se,
  log = TRUE,
  interval = TRUE,
  operator = c("-", "/"),
  level = 0.95,
  df = NULL
)
```

## Arguments

- est:

  a numeric vector of length 2 containing the two estimated ED values.

- se:

  a numeric vector of length 2 containing the two standard errors.

- log:

  logical indicating whether or not estimates and standard errors are on
  log scale.

- interval:

  logical indicating whether or not a confidence interval should be
  returned.

- operator:

  character string taking one of the two values `"-"` (default) or `"/"`
  corresponding to a comparison based on the difference or the ratio.

- level:

  numeric value giving the confidence level.

- df:

  numeric value specifying the degrees of freedom for the percentile
  used in the confidence interval (optional). By default confidence
  interval relies on the normal distribution.

## Value

A matrix with the estimated difference or ratio and the associated
standard error and the resulting confidence interval (unless not
requested).

## Author

Christian Ritz

## Examples

``` r
## Comparing ED50 values as a ratio
comped(c(28.396, 65.573), c(1.875, 5.619), log = FALSE, operator = "/")
#> 
#> Estimated ratio of effective doses
#> 
#>      Estimate Std. Error    Lower  Upper
#> [1,] 0.433044   0.046847 0.341226 0.5249

## Comparing ED50 values as a difference
comped(c(28.396, 65.573), c(1.875, 5.619), log = FALSE, operator = "-")
#> 
#> Estimated difference of effective doses
#> 
#>      Estimate Std. Error    Lower   Upper
#> [1,] -37.1770     5.9236 -48.7870 -25.567
```
