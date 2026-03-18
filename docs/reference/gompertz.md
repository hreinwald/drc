# Gompertz dose-response or growth curve model

Provides a very general way of specifying the mean function of the
decreasing or increasing Gompertz dose-response or growth curve models.

## Usage

``` r
gompertz(
  fixed = c(NA, NA, NA, NA),
  names = c("b", "c", "d", "e"),
  method = c("1", "2", "3", "4"),
  ssfct = NULL,
  fctName,
  fctText
)
```

## Arguments

- fixed:

  numeric vector. Specifies which parameters are fixed and at what value
  they are fixed. NAs for parameters that are not fixed.

- names:

  vector of character strings giving the names of the parameters (should
  not contain ":"). The order of the parameters is: b, c, d, e.

- method:

  character string indicating the self starter function to use.

- ssfct:

  a self starter function to be used.

- fctName:

  optional character string used internally by convenience functions.

- fctText:

  optional character string used internally by convenience functions.

## Value

A list containing the non-linear function, the self starter function and
the parameter names.

## Details

The Gompertz model is given by the mean function \$\$f(x) = c +
(d-c)(\exp(-\exp(b(x-e))))\$\$

If \\b\<0\\ the mean function is increasing; it is decreasing for
\\b\>0\\.

## References

Seber, G. A. F. and Wild, C. J. (1989) *Nonlinear Regression*, New York:
Wiley & Sons (p. 331).

## See also

The Weibull model
[`weibull2`](https://hreinwald.github.io/drc/reference/weibull2.md) is
closely related to the Gompertz model.

## Author

Christian Ritz
