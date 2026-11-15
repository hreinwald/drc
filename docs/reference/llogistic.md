# The log-logistic function

A very general way of specifying log-logistic models under various
constraints on parameters.

## Usage

``` r
llogistic(
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

  numeric vector of length 5, specifying fixed parameters (use NA for
  non-fixed parameters).

- names:

  character vector of length 5, specifying the names of the parameters:
  b, c, d, e, f.

- method:

  character string indicating the self starter function to use.

- ssfct:

  a self starter function to be used.

- fctName:

  optional character string used internally.

- fctText:

  optional character string used internally.

## Value

A list containing the nonlinear function, the self starter function, and
the parameter names.

## Details

The five-parameter log-logistic function is given by the expression
\$\$f(x) = c + \frac{d-c}{(1+\exp(b(\log(x)-\log(e))))^f}\$\$

The parameters are: `b`, the slope around `e` (positive for decreasing
curves, negative for increasing curves); `c`, the lower asymptote; `d`,
the upper asymptote (the response of the untreated control for a
decreasing curve); `e`, the ED50, i.e. the dose half-way between `c` and
`d` when `f = 1`; and `f`, the asymmetry parameter (`f = 1` gives the
symmetric model). Fixing parameters via `fixed` yields the convenience
models [`LL.2`](https://hreinwald.github.io/drc/reference/LL.2.md),
[`LL.3`](https://hreinwald.github.io/drc/reference/LL.3.md),
[`LL.3u`](https://hreinwald.github.io/drc/reference/LL.3u.md),
[`LL.4`](https://hreinwald.github.io/drc/reference/LL.4.md) and
[`LL.5`](https://hreinwald.github.io/drc/reference/LL.5.md).

## References

Finney, D. J. (1979).

Seber, G. A. F. and Wild, C. J. (1989).

## See also

[`LL.2`](https://hreinwald.github.io/drc/reference/LL.2.md),
[`LL.3`](https://hreinwald.github.io/drc/reference/LL.3.md),
[`LL.4`](https://hreinwald.github.io/drc/reference/LL.4.md),
[`LL.5`](https://hreinwald.github.io/drc/reference/LL.5.md)

## Author

Christian Ritz

## Examples

``` r
## Four-parameter log-logistic model (f fixed at 1): fit with drm() and plot
m1 <- drm(rootl ~ conc, data = ryegrass, fct = llogistic(fixed = c(NA, NA, NA, NA, 1)))
plot(m1, type = "all", main = "Four-parameter log-logistic model")


## Effect of the lower asymptote c: c free, fixed at 0 or fixed at 3
m.c0 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(fixed = c(NA, 0, NA, NA)))
m.c3 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(fixed = c(NA, 3, NA, NA)))
plot(m1, type = "all", main = "Lower asymptote c")
plot(m.c0, type = "none", add = TRUE, col = "darkgreen")
plot(m.c3, type = "none", add = TRUE, col = "red")
legend("topright", c("c free", "c = 0", "c = 3"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")


## Effect of the upper asymptote d: d free, fixed at 6 or fixed at 9
m.d6 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(fixed = c(NA, NA, 6, NA)))
m.d9 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(fixed = c(NA, NA, 9, NA)))
plot(m1, type = "all", main = "Upper asymptote d")
plot(m.d6, type = "none", add = TRUE, col = "darkgreen")
plot(m.d9, type = "none", add = TRUE, col = "red")
legend("topright", c("d free", "d = 6", "d = 9"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")

```
