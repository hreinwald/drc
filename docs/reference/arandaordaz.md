# Asymptotic Regression Model

The base function for the asymptotic regression model, providing the
mean function and self starter for a three-parameter model.

## Usage

``` r
arandaordaz(fixed = c(NA, NA, NA), names = c("a", "b", "c"), fctName, fctText)
```

## Arguments

- fixed:

  numeric vector. Specifies which parameters are fixed and at what value
  they are fixed. Use `NA` for parameters that are not fixed. Must be of
  length 3 (order: lower limit `a`, upper limit `b`, rate `c`).

- names:

  character vector of length 3 giving the names of the parameters
  (should not contain ":").

- fctName:

  optional character string used internally by convenience functions.
  Defaults to `"arandaordaz"` if not provided.

- fctText:

  optional character string used internally by convenience functions.
  Defaults to `"Asymptotic regression"` if not provided.

## Value

A list of class `drcMean` with the following components:

- fct:

  The mean function taking arguments `dose` and `parm`.

- ssfct:

  Self-starter function for generating initial parameter estimates from
  data.

- names:

  Character vector of non-fixed parameter names.

- deriv1:

  Reserved first derivative slot (currently `NULL`).

- deriv2:

  Reserved second derivative slot (currently `NULL`).

- derivx:

  Reserved derivative-with-respect-to-x slot (currently `NULL`).

- edfct:

  Function for calculating effective dose (ED) values and their
  derivatives.

- inversion:

  Inverse mean function for back-calculating dose from response.

- name:

  Character string identifying the model function name.

- text:

  Character string with a human-readable model description.

- noParm:

  Integer giving the number of non-fixed parameters.

## Details

The asymptotic regression model is a three-parameter model with mean
function:

\$\$f(x) = a + (b-a)(1-\exp(-c x))\$\$

The parameter \\a\\ is the lower limit (the response at \\x=0\\), \\b\\
is the upper limit (asymptote) and \\c\>0\\ is the rate of the increase.
The default `names` are `c("a", "b", "c")`. The equivalent models
[`AR.3`](https://hreinwald.github.io/drc/reference/AR.3.md) and
[`AR.2`](https://hreinwald.github.io/drc/reference/AR.2.md)
(parameterised by the scale \\e = 1/c\\) are usually preferred.

## See also

[`AR.2`](https://hreinwald.github.io/drc/reference/AR.2.md),
[`AR.3`](https://hreinwald.github.io/drc/reference/AR.3.md),
[`EXD.2`](https://hreinwald.github.io/drc/reference/EXD.2.md),
[`EXD.3`](https://hreinwald.github.io/drc/reference/EXD.3.md)

## Author

Christian Ritz, Hannes Reinwald

## Examples

``` r
## Simulated increasing data: fit with drm() and plot
dose <- rep(c(0.5, 1, 2, 4, 8, 16, 32), each = 3)
mu <- arandaordaz()$fct(dose, matrix(c(1, 8, 0.15), length(dose), 3, byrow = TRUE))
dat <- data.frame(dose = dose, response = mu + 0.15 * sin(seq_along(dose) * 2.3))
m1 <- drm(response ~ dose, data = dat, fct = arandaordaz())
plot(m1, type = "all", main = "Asymptotic regression")


## Effect of the lower limit a (fixed at 0) and the upper asymptote b (fixed at 6)
m.a0 <- drm(response ~ dose, data = dat, fct = arandaordaz(fixed = c(0, NA, NA)))
m.b6 <- drm(response ~ dose, data = dat, fct = arandaordaz(fixed = c(NA, 6, NA)))
plot(m1, type = "all", main = "Fixed lower limit and asymptote")
plot(m.a0, type = "none", add = TRUE, col = "darkgreen")
plot(m.b6, type = "none", add = TRUE, col = "red")
legend("bottomright", c("free", "a = 0", "b = 6"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")

```
