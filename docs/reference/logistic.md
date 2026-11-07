# The general asymmetric five-parameter logistic model

The five-parameter logistic model given by the expression \$\$f(x) = c +
\frac{d - c}{(1 + \exp(b(x - e)))^f}\$\$

## Usage

``` r
logistic(
  fixed = c(NA, NA, NA, NA, NA),
  names = c("b", "c", "d", "e", "f"),
  method = c("1", "2", "3", "4"),
  ssfct = NULL,
  fctName,
  fctText
)
```

## Arguments

- fixed:

  numeric vector of length 5. Specifies which parameters are fixed and
  at what value they are fixed. `NA` indicates that the corresponding
  parameter is not fixed.

- names:

  character vector of length 5 giving the names of the parameters
  `(b, c, d, e, f)`. Default is `c("b", "c", "d", "e", "f")`.

- method:

  character string indicating the self starter function to use (`"1"`,
  `"2"`, `"3"`, or `"4"`).

- ssfct:

  a self starter function to be used. If `NULL` (default), a built-in
  self starter is selected via `method`.

- fctName:

  optional character string used internally to overwrite the function
  name.

- fctText:

  optional character string used internally to overwrite the description
  text.

## Value

A list of class `"Boltzmann"` containing the nonlinear function, self
starter function, and parameter names.

## Details

This model differs from the log-logistic in that it uses `x` directly
rather than `log(x)`. It is sometimes referred to as the Boltzmann
model.

The parameters are: `b`, the slope (positive for a decreasing curve,
negative for an increasing curve); `c`, the lower asymptote; `d`, the
upper asymptote; `e`, the dose at the inflection point (the ED50 when
`f = 1`); and `f`, the asymmetry parameter. Effective doses from
[`ED`](https://hreinwald.github.io/drc/reference/ED.md) follow the same
convention as for the log-logistic models: for a decreasing curve
ED\\p\\ is the dose giving a response \\p\\\\ down to `c`; for an
increasing curve it is \\p\\\\ Absolute response levels
(`type = "absolute"`) are supported as well.

## See also

[`L.3`](https://hreinwald.github.io/drc/reference/L.3.md),
[`L.4`](https://hreinwald.github.io/drc/reference/L.4.md),
[`L.5`](https://hreinwald.github.io/drc/reference/L.5.md),
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md)

## Author

Christian Ritz

## Examples

``` r
## Fit with drm() and plot
m1 <- drm(rootl ~ conc, data = ryegrass, fct = logistic(fixed = c(NA, NA, NA, NA, 1)))
plot(m1, type = "all", main = "Four-parameter logistic model")


## Effect of the lower asymptote c: c free, fixed at 0 or fixed at 3
m.c0 <- drm(rootl ~ conc, data = ryegrass, fct = L.4(fixed = c(NA, 0, NA, NA)))
m.c3 <- drm(rootl ~ conc, data = ryegrass, fct = L.4(fixed = c(NA, 3, NA, NA)))
plot(m1, type = "all", main = "Lower asymptote c")
plot(m.c0, type = "none", add = TRUE, col = "darkgreen")
plot(m.c3, type = "none", add = TRUE, col = "red")
legend("topright", c("c free", "c = 0", "c = 3"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")


## Effect of the upper asymptote d: d free, fixed at 6 or fixed at 9
m.d6 <- drm(rootl ~ conc, data = ryegrass, fct = L.4(fixed = c(NA, NA, 6, NA)))
m.d9 <- drm(rootl ~ conc, data = ryegrass, fct = L.4(fixed = c(NA, NA, 9, NA)))
plot(m1, type = "all", main = "Upper asymptote d")
plot(m.d6, type = "none", add = TRUE, col = "darkgreen")
plot(m.d9, type = "none", add = TRUE, col = "red")
legend("topright", c("d free", "d = 6", "d = 9"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")


## Effective doses
ED(m1, c(10, 50, 90))
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error
#> e:10  1.08166    0.51209
#> e:50  3.06924    0.19638
#> e:90  5.05682    0.38615
```
