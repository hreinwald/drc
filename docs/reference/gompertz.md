# Gompertz dose-response or growth curve model

Provides a very general way of specifying the mean function of the
decreasing or increasing Gompertz dose-response or growth curve models.

## Usage

``` r
gompertz(
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

  numeric vector. Specifies which parameters are fixed and at what value
  they are fixed. NAs for parameters that are not fixed.

- names:

  vector of character strings giving the names of the parameters (should
  not contain ":"). The order of the parameters is: b, c, d, e.

- method:

  character string indicating the self starter function to use.

- ssfct:

  a self starter function to be used.

- fctName:

  optional character string used internally by convenience functions.

- fctText:

  optional character string used internally by convenience functions.

## Value

A list containing the non-linear function, the self starter function and
the parameter names.

## Details

The Gompertz model is given by the mean function \$\$f(x) = c +
(d-c)(\exp(-\exp(b(x-e))))\$\$

If \\b\<0\\ the mean function is increasing; it is decreasing for
\\b\>0\\.

Note that the model is a function of the dose \\x\\ itself (not of
\\\log(x)\\). The parameters are: `b`, the slope; `c`, the lower
asymptote; `d`, the upper asymptote; and `e`, the location parameter
(the dose at which the response is \\c + (d-c)/\exp(1)\\, so it is not
the ED50). Effective doses are obtained with
[`ED`](https://hreinwald.github.io/drc/reference/ED.md): \\ED_p = e +
\log(-\log(1-p/100))/b\\ for a decreasing curve, with standard errors
from the delta method.

## References

Seber, G. A. F. and Wild, C. J. (1989) *Nonlinear Regression*, New York:
Wiley & Sons (p. 331).

## See also

The Weibull model
[`weibull2`](https://hreinwald.github.io/drc/reference/weibull2.md) is
closely related to the Gompertz model.

## Author

Christian Ritz

## Examples

``` r
## Fit with drm() and plot
m1 <- drm(rootl ~ conc, data = ryegrass, fct = gompertz())
plot(m1, type = "all", main = "Four-parameter Gompertz model")


## Effect of the lower asymptote c: c free, fixed at 0 or fixed at 3
m.c0 <- drm(rootl ~ conc, data = ryegrass, fct = G.4(fixed = c(NA, 0, NA, NA)))
m.c3 <- drm(rootl ~ conc, data = ryegrass, fct = G.4(fixed = c(NA, 3, NA, NA)))
plot(m1, type = "all", main = "Lower asymptote c")
plot(m.c0, type = "none", add = TRUE, col = "darkgreen")
plot(m.c3, type = "none", add = TRUE, col = "red")
legend("topright", c("c free", "c = 0", "c = 3"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")


## Effect of the upper asymptote d: d free, fixed at 6 or fixed at 9
m.d6 <- drm(rootl ~ conc, data = ryegrass, fct = G.4(fixed = c(NA, NA, 6, NA)))
m.d9 <- drm(rootl ~ conc, data = ryegrass, fct = G.4(fixed = c(NA, NA, 9, NA)))
plot(m1, type = "all", main = "Upper asymptote d")
plot(m.d6, type = "none", add = TRUE, col = "darkgreen")
plot(m.d9, type = "none", add = TRUE, col = "red")
legend("topright", c("d free", "d = 6", "d = 9"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")


## Effective doses with standard errors
ED(m1, c(10, 50, 90))
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error
#> e:10  0.85872    0.72133
#> e:50  3.11669    0.22075
#> e:90  4.55565    0.25033
```
