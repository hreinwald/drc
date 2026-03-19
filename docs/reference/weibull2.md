# The four-parameter Weibull (type 2) model

Provides a general framework for the four-parameter Weibull type 2 model
given by the equation \$\$f(x) = c + (d - c)(1 - \exp(-\exp(b(\log(x) -
\log(e)))))\$\$

## Usage

``` r
weibull2(
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

  numeric vector of length 4, specifying fixed parameters (use `NA` for
  parameters that should be estimated).

- names:

  character vector of length 4 giving the names of the parameters
  (default `c("b", "c", "d", "e")`).

- method:

  character string indicating the self starter method to use for
  obtaining starting values. One of `"1"` (default), `"2"`, `"3"`, or
  `"4"`. See Details.

- ssfct:

  a self starter function. If `NULL` (default), a built-in self starter
  is used based on `method`.

- fctName:

  optional character string used internally for the function name.

- fctText:

  optional character string used internally for the function
  description.

## Value

A list containing the nonlinear function, self starter function, and
parameter names. The list has class `"Weibull-2"`.

## Details

The `method` argument determines how starting values for the parameters
`b` and `e` are estimated (the starting values for `c` and `d` are
always based on the range of the response values). Four methods are
available:

- `"1"` (default):

  Linear regression on transformed data. Applies a complementary log-log
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

[`weibull1`](https://hreinwald.github.io/drc/reference/weibull1.md),
[`W2.2`](https://hreinwald.github.io/drc/reference/W2.2.md),
[`W2.3`](https://hreinwald.github.io/drc/reference/W2.3.md),
[`W2.4`](https://hreinwald.github.io/drc/reference/W2.4.md)

## Author

Christian Ritz
