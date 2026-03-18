# Three-Phase Dose-Response Model

A ten-parameter dose-response model combining three log-logistic
components, extending the two-phase model
([`twophase`](https://hreinwald.github.io/drc/reference/twophase.md))
for describing even more complex dose-response patterns.

## Usage

``` r
threephase(
  fixed = c(NA, NA, NA, NA, NA, NA, NA, NA, NA, NA),
  names = c("b1", "c1", "d1", "e1", "b2", "d2", "e2", "b3", "d3", "e3"),
  fctName,
  fctText
)
```

## Arguments

- fixed:

  numeric vector specifying which parameters are fixed and at what value
  they are fixed. NAs are used for parameters that are not fixed.

- names:

  a vector of character strings giving the names of the parameters
  (should not contain ":"). The default is reasonable.

- fctName:

  optional character string used internally by convenience functions.

- fctText:

  optional character string used internally by convenience functions.

## Value

A list containing the nonlinear function, the self starter function, and
the parameter names.

## Details

The model function is the sum of a four-parameter log-logistic model and
two three-parameter log-logistic models:

\$\$f(x) = \mathrm{LL.4}(x; b1, c1, d1, e1) + \mathrm{LL.3}(x; b2, d2,
e2) + \mathrm{LL.3}(x; b3, d3, e3)\$\$

## See also

[`twophase`](https://hreinwald.github.io/drc/reference/twophase.md),
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md)

## Author

Christian Ritz
