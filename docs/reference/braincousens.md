# The Brain-Cousens hormesis models

`braincousens` provides a very general way of specifying Brain-Cousens'
modified log-logistic model for describing hormesis, under various
constraints on the parameters.

## Usage

``` r
braincousens(
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

  numeric vector. Specifies which parameters are fixed and at what value
  they are fixed. NAs for parameters that are not fixed.

- names:

  a vector of character strings giving the names of the parameters
  (should not contain ":"). The order of the parameters is: b, c, d, e,
  f.

- method:

  character string indicating the self starter function to use.

- ssfct:

  a self starter function to be used.

- fctName:

  optional character string used internally by convenience functions.

- fctText:

  optional character string used internally by convenience functions.

## Value

A list containing the non-linear function, the self starter function,
the parameter names and additional model specific objects.

## Details

The Brain-Cousens model is given by the expression \$\$f(x) = c +
\frac{d-c+fx}{1+\exp(b(\log(x)-\log(e)))}\$\$ which is a five-parameter
model.

The parameters are: `b`, the slope; `c`, the lower asymptote; `d`, the
upper asymptote (the response of the untreated control is `d`); `e`, the
location of the decline; and `f`, the hormesis parameter (the initial
slope of the stimulation at low doses; `f = 0` gives the log-logistic
model [`LL.4`](https://hreinwald.github.io/drc/reference/LL.4.md)).
Effective doses are computed relative to the asymptotes `c` and `d` (see
[`ED`](https://hreinwald.github.io/drc/reference/ED.md), use
`bound = FALSE` if levels beyond the range are needed) and the maximal
response with [`MAX`](https://hreinwald.github.io/drc/reference/MAX.md).

## References

Brain, P. and Cousens, R. (1989) An equation to describe dose responses
where there is stimulation of growth at low doses, *Weed Research*,
**29**, 93–96.

## See also

[`BC.4`](https://hreinwald.github.io/drc/reference/BC.4.md),
[`BC.5`](https://hreinwald.github.io/drc/reference/BC.5.md),
[`drm`](https://hreinwald.github.io/drc/reference/drm.md)

## Author

Christian Ritz

## Examples

``` r
## Fit with drm() and plot
m1 <- drm(weight ~ conc, data = lettuce, fct = braincousens())
plot(m1, type = "all", main = "Brain-Cousens hormesis model")

MAX(m1)
#>     Dose Response
#> 1 0.7203   1.2008

## Effect of the lower asymptote c: c free, fixed at 0 or fixed at 0.3
m.c0 <- drm(weight ~ conc, data = lettuce, fct = BC.5(fixed = c(NA, 0, NA, NA, NA)))
m.c3 <- drm(weight ~ conc, data = lettuce, fct = BC.5(fixed = c(NA, 0.3, NA, NA, NA)))
plot(m1, type = "all", main = "Lower asymptote c")
plot(m.c0, type = "none", add = TRUE, col = "darkgreen")
plot(m.c3, type = "none", add = TRUE, col = "red")
legend("topright", c("c free", "c = 0", "c = 0.3"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")


## Effect of the upper asymptote d: d free, fixed at 0.8 or fixed at 1.2
m.d8 <- drm(weight ~ conc, data = lettuce, fct = BC.5(fixed = c(NA, NA, 0.8, NA, NA)))
m.d12 <- drm(weight ~ conc, data = lettuce, fct = BC.5(fixed = c(NA, NA, 1.2, NA, NA)))
plot(m1, type = "all", main = "Upper asymptote d")
plot(m.d8, type = "none", add = TRUE, col = "darkgreen")
plot(m.d12, type = "none", add = TRUE, col = "red")
legend("topright", c("d free", "d = 0.8", "d = 1.2"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")

```
