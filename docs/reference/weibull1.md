# The four-parameter Weibull type 1 model

The general Weibull type 1 model for fitting dose-response data.

## Usage

``` r
weibull1(
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

  numeric vector of length 4. Specifies which parameters are fixed and
  at what value. Use `NA` for parameters that are not fixed.

- names:

  character vector of length 4 giving the names of the parameters `b`,
  `c`, `d`, and `e`.

- method:

  character string indicating the self starter function to use (`"1"`,
  `"2"`, `"3"`, or `"4"`).

- ssfct:

  a self starter function to be used. If `NULL` (default), the built-in
  self starter is used.

- fctName:

  optional character string used internally for the function name.

- fctText:

  optional character string used internally for the function text
  description.

## Value

A list of class `Weibull-1` containing the nonlinear function, self
starter function, and parameter names.

## Details

The four-parameter Weibull type 1 model is given by the expression
\$\$f(x) = c + (d - c) \exp(-\exp(b(\log(x) - \log(e))))\$\$

The model is sometimes also called the Gompertz model.

## References

Seber, G. A. F. and Wild, C. J. (1989) *Nonlinear Regression*, New York:
Wiley & Sons (pp. 338–339).

## See also

[`W1.2`](https://hreinwald.github.io/drc/reference/W1.2.md),
[`W1.3`](https://hreinwald.github.io/drc/reference/W1.3.md),
[`W1.4`](https://hreinwald.github.io/drc/reference/W1.4.md),
[`weibull2`](https://hreinwald.github.io/drc/reference/weibull2.md)

## Author

Christian Ritz
