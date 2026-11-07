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

The parameters are: `b`, the slope (positive for a decreasing curve);
`c`, the lower asymptote; `d`, the upper asymptote (the control response
of a decreasing curve); and `e`, the dose at which the curve has dropped
to `c + (d - c) / exp(1)`, about 37\\ way from `d` to `c`. Note that `e`
is not the ED50; use
[`ED`](https://hreinwald.github.io/drc/reference/ED.md) for effective
doses. The curve is asymmetric on the log-dose axis.

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

## Examples

``` r
## Fit with drm() and plot
m1 <- drm(rootl ~ conc, data = ryegrass, fct = weibull1())
plot(m1, type = "all", main = "Four-parameter Weibull type 1 model")


## Effect of the lower asymptote c: c free, fixed at 0 or fixed at 3
m.c0 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4(fixed = c(NA, 0, NA, NA)))
m.c3 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4(fixed = c(NA, 3, NA, NA)))
plot(m1, type = "all", main = "Lower asymptote c")
plot(m.c0, type = "none", add = TRUE, col = "darkgreen")
plot(m.c3, type = "none", add = TRUE, col = "red")
legend("topright", c("c free", "c = 0", "c = 3"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")


## Effect of the upper asymptote d: d free, fixed at 6 or fixed at 9
m.d6 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4(fixed = c(NA, NA, 6, NA)))
m.d9 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4(fixed = c(NA, NA, 9, NA)))
plot(m1, type = "all", main = "Upper asymptote d")
plot(m.d6, type = "none", add = TRUE, col = "darkgreen")
plot(m.d9, type = "none", add = TRUE, col = "red")
legend("topright", c("d free", "d = 6", "d = 9"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")

```
