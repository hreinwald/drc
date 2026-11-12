# Fractional polynomial-logistic dose-response model

Model function for specifying dose-response models that are a
combination of a logistic model and an appropriate class of fractional
polynomials.

## Usage

``` r
fplogistic(
  p1,
  p2,
  fixed = c(NA, NA, NA, NA),
  names = c("b", "c", "d", "e"),
  method = c("1", "2", "3", "4"),
  ssfct = NULL,
  fctName,
  fctText
)
```

## Arguments

- p1:

  numeric denoting the negative power of log(dose+1) in the fractional
  polynomial.

- p2:

  numeric denoting the positive power of log(dose+1) in the fractional
  polynomial.

- fixed:

  numeric vector. Specifies which parameters are fixed and at what value
  they are fixed. NAs for parameters that are not fixed.

- names:

  a vector of character strings giving the names of the parameters
  (should not contain ":"). The order of the parameters is: b, c, d, e.

- method:

  character string indicating the self starter function to use.

- ssfct:

  a self starter function to be used.

- fctName:

  optional character string used internally by convenience functions.

- fctText:

  optional character string used internally by convenience functions.

## Value

A list containing the nonlinear function, the self starter function and
the parameter names.

## Details

The fractional polynomial dose-response models introduced by Namata et
al. (2008) are implemented using the logistic model as base: \$\$f(x) =
c + \frac{d-c}{1+\exp(b \log(x+1)^{p1} + e \log(x+1)^{p2})}\$\$ where
`c` is the lower asymptote, `d` the upper asymptote and `b` and `e` are
the coefficients of the two fractional polynomial terms (`e` is not an
ED50). The powers `p1` and `p2` are fixed by the user (typically from
\\\\-2, -1, -0.5, 0, 0.5, 1, 2, 3\\\\).

## References

Namata, Harriet and Aerts, Marc and Faes, Christel and Teunis, Peter
(2008) Model Averaging in Microbial Risk Assessment Using Fractional
Polynomials, *Risk Analysis* **28**, 891–905.

## See also

[`FPL.4`](https://hreinwald.github.io/drc/reference/FPL.4.md),
[`maED`](https://hreinwald.github.io/drc/reference/maED.md),
[`drm`](https://hreinwald.github.io/drc/reference/drm.md)

## Author

Christian Ritz

## Examples

``` r
## Fit with drm() and plot
m1 <- drm(rootl ~ conc, data = ryegrass, fct = fplogistic(p1 = -1, p2 = 1))
plot(m1, type = "all", main = "Fractional polynomial-logistic model")


## Effect of the lower asymptote c and the upper asymptote d
m.c0 <- drm(rootl ~ conc, data = ryegrass, fct = FPL.4(-1, 1, fixed = c(NA, 0, NA, NA)))
m.d6 <- drm(rootl ~ conc, data = ryegrass, fct = FPL.4(-1, 1, fixed = c(NA, NA, 6, NA)))
plot(m1, type = "all", main = "Fixed asymptotes")
plot(m.c0, type = "none", add = TRUE, col = "darkgreen")
plot(m.d6, type = "none", add = TRUE, col = "red")
legend("topright", c("free", "c = 0", "d = 6"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")

```
