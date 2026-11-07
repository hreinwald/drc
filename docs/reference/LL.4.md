# Four-parameter log-logistic function

A four-parameter log-logistic function.

## Usage

``` r
LL.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)

l4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 4, specifying fixed parameters (use NA for
  non-fixed parameters).

- names:

  character vector of length 4, specifying the names of the parameters
  (default: b, c, d, e).

- ...:

  additional arguments to
  [`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Value

See
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Details

The four-parameter log-logistic function is given by the expression
\$\$f(x) = c + \frac{d-c}{1+\exp(b(\log(x)-\log(e)))}\$\$

## See also

[`LL.3`](https://hreinwald.github.io/drc/reference/LL.3.md),
[`LL.5`](https://hreinwald.github.io/drc/reference/LL.5.md),
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md)

## Author

Christian Ritz and Jens C. Streibig

## Examples

``` r
## Several curves fitted with drm() and plotted
spinach.m1 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4())
plot(spinach.m1, main = "LL.4, one curve per CURVE level")


## Effect of the asymptotes on a single curve (ryegrass): fix c at 0 / d at 6
m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
m.c0 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(fixed = c(NA, 0, NA, NA)))
m.d6 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(fixed = c(NA, NA, 6, NA)))
plot(m1, type = "all", main = "LL.4: free and fixed asymptotes")
plot(m.c0, type = "none", add = TRUE, col = "darkgreen")
plot(m.d6, type = "none", add = TRUE, col = "red")
legend("topright", c("free", "c = 0", "d = 6"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")

```
