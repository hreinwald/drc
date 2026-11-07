# Gamma Dose-Response Model

A four-parameter dose-response model derived from the cumulative
distribution function of the gamma distribution. Only suitable for
increasing dose-response data.

## Usage

``` r
gammadr(
  fixed = c(NA, NA, NA, NA),
  names = c("b", "c", "d", "e"),
  fctName,
  fctText
)
```

## Arguments

- fixed:

  numeric vector specifying which parameters are fixed and at what value
  they are fixed. NAs are used for parameters that are not fixed.

- names:

  a vector of character strings giving the names of the parameters
  (should not contain ":"). The default is reasonable.

- fctName:

  optional character string used internally by convenience functions.

- fctText:

  optional character string used internally by convenience functions.

## Value

A list containing the nonlinear function, the self starter function, and
the parameter names.

## Details

Following Wheeler and Bailer (2009) the model function is:

\$\$f(x) = c + (d-c) \cdot \mathrm{pgamma}(b \cdot x, e, 1)\$\$

The parameters are: `b`, the rate (scale) of the dose; `c`, the lower
asymptote (the response at dose 0); `d`, the upper asymptote (approached
at high doses); and `e`, the shape of the gamma distribution. Effective
doses are not available for this model;
[`predict`](https://rdrr.io/r/stats/predict.html) uses the analytic
derivatives for standard errors.

## References

Wheeler, M. W., Bailer, A. J. (2009) Comparing model averaging with
other model selection strategies for benchmark dose estimation,
*Environmental and Ecological Statistics*, **16**, 37–51.

## See also

[`drm`](https://hreinwald.github.io/drc/reference/drm.md),
[`multi2`](https://hreinwald.github.io/drc/reference/multi2.md)

## Author

Christian Ritz

## Examples

``` r
## Simulated increasing dose-response data: fit with drm() and plot
dose <- rep(c(0.5, 1, 2, 4, 8, 16, 32), each = 3)
mu <- gammadr()$fct(dose, matrix(c(0.2, 1, 9, 2.5), length(dose), 4, byrow = TRUE))
dat <- data.frame(dose = dose, response = mu + 0.25 * sin(seq_along(dose) * 2.3))
m1 <- drm(response ~ dose, data = dat, fct = gammadr())
plot(m1, type = "all", main = "Gamma dose-response model")


## Effect of the lower asymptote c and the upper asymptote d
m.c0 <- drm(response ~ dose, data = dat, fct = gammadr(fixed = c(NA, 0, NA, NA)))
m.d7 <- drm(response ~ dose, data = dat, fct = gammadr(fixed = c(NA, NA, 7, NA)))
plot(m1, type = "all", main = "Fixed asymptotes")
plot(m.c0, type = "none", add = TRUE, col = "darkgreen")
plot(m.d7, type = "none", add = TRUE, col = "red")
legend("bottomright", c("free", "c = 0", "d = 7"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")

```
