# Gamma Dose-Response Model

A four-parameter dose-response model derived from the cumulative
distribution function of the gamma distribution. Only suitable for
increasing dose-response data.

## Usage

``` r
gammadr(
  fixed = c(NA, NA, NA, NA),
  names = c("b", "c", "d", "e"),
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

Following Wheeler and Bailer (2009) the model function is:

\$\$f(x) = c + (d-c) \cdot \mathrm{pgamma}(b \cdot x, e, 1)\$\$

## References

Wheeler, M. W., Bailer, A. J. (2009) Comparing model averaging with
other model selection strategies for benchmark dose estimation,
*Environmental and Ecological Statistics*, **16**, 37–51.

## Author

Christian Ritz
