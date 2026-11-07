# Five-Parameter Log-Logistic Model with log(ED50) as Parameter

A five-parameter log-logistic model where the ED50 is parameterised on
the log scale. The mean function is: \$\$f(x) = c + \frac{d - c}{(1 +
\exp(b(\log(x) - e)))^f}\$\$ where `e` is the logarithm of the ED50 (not
exponentiated).

## Usage

``` r
llogistic2(
  fixed = c(NA, NA, NA, NA, NA),
  names = c("b", "c", "d", "e", "f"),
  ss = c("1", "2", "3"),
  ssfct = NULL,
  fctName,
  fctText
)
```

## Arguments

- fixed:

  numeric vector of length 5. Specifies which parameters are fixed and
  at what value. Use `NA` for parameters that should be estimated.

- names:

  character vector of length 5 giving the names of the parameters `b`,
  `c`, `d`, `e`, and `f`.

- ss:

  character string indicating the self-starter version to use. One of
  `"1"` (default), `"2"`, or `"3"`.

- ssfct:

  optional self-starter function. If provided, overrides the built-in
  self-starter selected by `ss`.

- fctName:

  optional character string specifying the name of the function.

- fctText:

  optional character string providing a short description of the
  function.

## Value

A list of class `"llogistic"` containing the nonlinear function,
self-starter function, parameter names, and related helper functions.

## Details

Apart from the parameterisation of `e` the model is identical to
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md):
`b` is the slope, `c` the lower asymptote, `d` the upper asymptote and
`f` the asymmetry parameter.

Because `e` is on the log scale,
[`ED`](https://hreinwald.github.io/drc/reference/ED.md) returns
**log(ED)** for these models. Use
`ED(object, respLev, interval = "fls")` to obtain the ED values (and
confidence limits) on the original dose scale; `interval = "fls"` is
only valid for models of this kind. Model averaging with
[`maED`](https://hreinwald.github.io/drc/reference/maED.md) is not
available for these models (use the dose-scale equivalents
[`LL.2`](https://hreinwald.github.io/drc/reference/LL.2.md), ...,
[`LL.5`](https://hreinwald.github.io/drc/reference/LL.5.md) instead).

## See also

[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md),
[`LL2.2`](https://hreinwald.github.io/drc/reference/LL2.2.md),
[`LL2.3`](https://hreinwald.github.io/drc/reference/LL2.3.md),
[`LL2.4`](https://hreinwald.github.io/drc/reference/LL2.4.md),
[`LL2.5`](https://hreinwald.github.io/drc/reference/LL2.5.md)

## Author

Christian Ritz

## Examples

``` r
## Four-parameter model with log(ED50) as parameter: fit with drm() and plot
m1 <- drm(rootl ~ conc, data = ryegrass, fct = llogistic2(fixed = c(NA, NA, NA, NA, 1)))
plot(m1, type = "all", main = "Log-logistic model with log(ED50) as parameter")


## ED values on the dose scale (interval = "fls" back-transforms from the log scale)
ED(m1, c(10, 50, 90), interval = "fls")
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error    Lower    Upper
#> e:10 1.463710   0.127602 1.121651 1.910083
#> e:50 3.057955   0.060737 2.694064 3.470996
#> e:90 6.388620   0.132282 4.848085 8.418679

## Effect of the lower asymptote c and the upper asymptote d
m.c0 <- drm(rootl ~ conc, data = ryegrass, fct = LL2.4(fixed = c(NA, 0, NA, NA)))
m.c3 <- drm(rootl ~ conc, data = ryegrass, fct = LL2.4(fixed = c(NA, 3, NA, NA)))
plot(m1, type = "all", main = "Lower asymptote c")
plot(m.c0, type = "none", add = TRUE, col = "darkgreen")
plot(m.c3, type = "none", add = TRUE, col = "red")
legend("topright", c("c free", "c = 0", "c = 3"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")


m.d6 <- drm(rootl ~ conc, data = ryegrass, fct = LL2.4(fixed = c(NA, NA, 6, NA)))
m.d9 <- drm(rootl ~ conc, data = ryegrass, fct = LL2.4(fixed = c(NA, NA, 9, NA)))
plot(m1, type = "all", main = "Upper asymptote d")
plot(m.d6, type = "none", add = TRUE, col = "darkgreen")
plot(m.d9, type = "none", add = TRUE, col = "red")
legend("topright", c("d free", "d = 6", "d = 9"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")

```
