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

The parameters are: `b`, the slope (negative for a decreasing curve,
positive for an increasing curve); `c`, the lower asymptote; `d`, the
upper asymptote; and `e`, the dose at which the response has covered
about 63\\ ED50; use
[`ED`](https://hreinwald.github.io/drc/reference/ED.md)). The curve is
asymmetric on the log-dose axis and complementary to
[`weibull1`](https://hreinwald.github.io/drc/reference/weibull1.md).

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

## Examples

``` r
## Fit with drm() and plot
m1 <- drm(rootl ~ conc, data = ryegrass, fct = weibull2())
plot(m1, type = "all", main = "Four-parameter Weibull type 2 model")


## Effect of the lower asymptote c: c free, fixed at 0 or fixed at 3
m.c0 <- drm(rootl ~ conc, data = ryegrass, fct = W2.4(fixed = c(NA, 0, NA, NA)))
m.c3 <- drm(rootl ~ conc, data = ryegrass, fct = W2.4(fixed = c(NA, 3, NA, NA)))
plot(m1, type = "all", main = "Lower asymptote c")
plot(m.c0, type = "none", add = TRUE, col = "darkgreen")
plot(m.c3, type = "none", add = TRUE, col = "red")
legend("topright", c("c free", "c = 0", "c = 3"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")


## Effect of the upper asymptote d: d free, fixed at 6 or fixed at 9
m.d6 <- drm(rootl ~ conc, data = ryegrass, fct = W2.4(fixed = c(NA, NA, 6, NA)))
m.d9 <- drm(rootl ~ conc, data = ryegrass, fct = W2.4(fixed = c(NA, NA, 9, NA)))
plot(m1, type = "all", main = "Upper asymptote d")
plot(m.d6, type = "none", add = TRUE, col = "darkgreen")
plot(m.d9, type = "none", add = TRUE, col = "red")
legend("topright", c("d free", "d = 6", "d = 9"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")

```
