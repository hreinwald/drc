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

  character string indicating the self starter function to use for
  obtaining starting values (`"1"` (default), `"2"`, `"3"`, or `"4"`).
  See Details.

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

The `method` argument determines how starting values for the parameters
`b` and `e` are estimated (the starting values for `c` and `d` are
always based on the range of the response values). Four methods are
available:

- `"1"` (default):

  Linear regression on transformed data. Applies a log-log
  transformation to the response and a log transformation to the dose,
  then fits a linear regression to estimate starting values for `b` and
  `e`.

- `"2"`:

  Anke's procedure. Estimates `e` by finding the dose at which the
  response crosses the midpoint between `c` and `d`, then estimates `b`
  as the median of back-calculated values.

- `"3"`:

  Stepwise approach. Identifies where the mean response crosses the
  midpoint between `c` and `d` and uses the corresponding dose as the
  starting value for `e`. The starting value for `b` is based on the
  sign of the slope at that point.

- `"4"`:

  Normolle's procedure. Uses the mean of the dose range as an initial
  estimate for `e`, then estimates `b` and `e` using median-based
  back-calculations.

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
