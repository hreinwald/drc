# Two-Phase Dose-Response Model

A seven-parameter dose-response model combining two log-logistic
components, useful for describing more complex dose-response patterns.

## Usage

``` r
twophase(
  fixed = c(NA, NA, NA, NA, NA, NA, NA),
  names = c("b1", "c1", "d1", "e1", "b2", "d2", "e2"),
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

Following Groot *et al* (1996) the two-phase model function is:

\$\$f(x) = c + \frac{d1-c}{1+\exp(b1(\log(x)-\log(e1)))} +
\frac{d2}{1+\exp(b2(\log(x)-\log(e2)))}\$\$

For each of the two phases, the parameters have the same interpretation
as in the ordinary log-logistic model: `b1`, `b2` are the slopes and
`e1`, `e2` the locations of the two phases; `c` (here `c1`) is the lower
asymptote of the whole curve and the upper asymptote is `d1 + d2`.
Starting values from the self starter may need to be replaced by
user-supplied `start` values. No effective doses
([`ED`](https://hreinwald.github.io/drc/reference/ED.md)) are available.

## References

Groot, J. C. J., Cone, J. W., Williams, B. A., Debersaques, F. M. A.,
Lantinga, E. A. (1996) Multiphasic analysis of gas production kinetics
for in vitro fermentation of ruminant feeds, *Animal Feed Science
Technology*, **64**, 77–89.

## See also

The basic component in the two-phase model is the log-logistic model
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Author

Christian Ritz

## Examples

``` r
## Simulated two-phase data: fit with drm() and plot (starting values supplied)
dose <- rep(10^seq(-1, 4, length.out = 12), each = 2)
parm <- c(b1 = 2, c1 = 0.5, d1 = 4, e1 = 3, b2 = 1.5, d2 = 4, e2 = 400)
mu <- twophase()$fct(dose, matrix(parm, length(dose), 7, byrow = TRUE))
dat <- data.frame(dose = dose, response = mu + 0.04 * sin(seq_along(dose) * 1.7))
m1 <- drm(response ~ dose, data = dat, fct = twophase(), start = parm)
plot(m1, type = "all", log = "x", main = "Two-phase model")

```
