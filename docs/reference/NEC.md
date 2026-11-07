# No Effect Concentration (NEC) dose-response model

The NEC model is a dose-response model with a threshold below which the
response is assumed constant and equal to the control response. It has
been proposed as an alternative to both the classical NOEC and the
regression-based EC/ED approach.

## Usage

``` r
NEC(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), fctName, fctText)
```

## Arguments

- fixed:

  numeric vector specifying which parameters are fixed and at what value
  they are fixed. NAs are used for parameters that are not fixed.

- names:

  a vector of character strings giving the names of the parameters
  (should not contain ":"). The default is reasonable (see under
  'Usage').

- fctName:

  optional character string used internally by convenience functions.

- fctText:

  optional character string used internally by convenience functions.

## Value

A list containing the nonlinear function, the self starter function and
the parameter names.

## Details

The NEC model function proposed by Pires et al (2002) is: \$\$f(x) = c +
(d-c)\exp(-b(x-e)I(x-e))\$\$ where \\I(x-e)\\ is the indicator function
equal to 0 for \\x\<=e\\ and 1 for \\x\>e\\.

The parameters are: `b`, the rate of the exponential decline above the
threshold; `c`, the lower asymptote (the response approached at high
doses); `d`, the upper asymptote, i.e. the constant response up to the
threshold; and `e`, the threshold itself, the no-effect concentration
(NEC), which is estimated directly as a model parameter.

Effective doses are not available for these models
([`ED`](https://hreinwald.github.io/drc/reference/ED.md) stops with an
error). Standard errors and intervals from
[`predict`](https://rdrr.io/r/stats/predict.html) are based on numerical
derivatives.

## References

Pires, A. M., Branco, J. A., Picado, A., Mendonca, E. (2002) Models for
the estimation of a 'no effect concentration', *Environmetrics*, **13**,
15–27.

## See also

[`NEC.2`](https://hreinwald.github.io/drc/reference/NEC.2.md),
[`NEC.3`](https://hreinwald.github.io/drc/reference/NEC.3.md),
[`NEC.4`](https://hreinwald.github.io/drc/reference/NEC.4.md),
[`drm`](https://hreinwald.github.io/drc/reference/drm.md)

## Author

Christian Ritz

## Examples

``` r
## Fit with drm() and plot; the parameter e is the NEC
nec.m1 <- drm(rootl ~ conc, data = ryegrass, fct = NEC.4())
summary(nec.m1)
#> 
#> Model fitted: NEC (4 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) 0.447556   0.072119  6.2058 4.623e-06 ***
#> c:(Intercept) 0.498164   0.204584  2.4350   0.02438 *  
#> d:(Intercept) 7.723988   0.173236 44.5865 < 2.2e-16 ***
#> e:(Intercept) 1.428649   0.160386  8.9076 2.132e-08 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.5197084 (20 degrees of freedom)
plot(nec.m1, type = "all", main = "NEC model")


## Effect of the lower asymptote c: c free, fixed at 0 or fixed at 3
nec.c0 <- drm(rootl ~ conc, data = ryegrass, fct = NEC.4(fixed = c(NA, 0, NA, NA)))
nec.c3 <- drm(rootl ~ conc, data = ryegrass, fct = NEC.4(fixed = c(NA, 3, NA, NA)))
plot(nec.m1, type = "all", main = "Lower asymptote c")
plot(nec.c0, type = "none", add = TRUE, col = "darkgreen")
plot(nec.c3, type = "none", add = TRUE, col = "red")
legend("topright", c("c free", "c = 0", "c = 3"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")


## Effect of the upper asymptote d: d free, fixed at 6 or fixed at 9
nec.d6 <- drm(rootl ~ conc, data = ryegrass, fct = NEC.4(fixed = c(NA, NA, 6, NA)))
nec.d9 <- drm(rootl ~ conc, data = ryegrass, fct = NEC.4(fixed = c(NA, NA, 9, NA)))
plot(nec.m1, type = "all", main = "Upper asymptote d")
plot(nec.d6, type = "none", add = TRUE, col = "darkgreen")
plot(nec.d9, type = "none", add = TRUE, col = "red")
legend("topright", c("d free", "d = 6", "d = 9"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")

```
