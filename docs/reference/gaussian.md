# Normal (Gaussian) biphasic dose-response model

Model function for fitting symmetric or skewed bell-shaped/biphasic
dose-response patterns using the Gaussian (normal distribution) model.

## Usage

``` r
gaussian(
  fixed = c(NA, NA, NA, NA, NA),
  names = c("b", "c", "d", "e", "f"),
  method = c("1", "2", "3", "4"),
  ssfct = NULL,
  fctName,
  fctText,
  loge = FALSE
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

- loge:

  logical indicating whether or not e or log(e) should be a parameter in
  the model. By default e is a model parameter.

## Value

The value returned is a list containing the nonlinear function, the self
starter function and the parameter names.

## Details

The model function is \$\$f(x) = c +
(d-c)\exp\left(-\frac{1}{2}\left\|\frac{x-e}{b}\right\|^f\right)\$\$ The
parameters are: `b`, the width (scale) of the bell; `c`, the lower
asymptote (the baseline response far from the peak); `d`, the response
at the peak (reached at the dose `e`, the location of the peak, or its
logarithm if `loge = TRUE`); and `f`, the shape (`f = 2` gives the
normal density). If `d` exceeds `c` the response peaks, otherwise it
dips. No effective doses
([`ED`](https://hreinwald.github.io/drc/reference/ED.md)) are available
for this model.

## See also

[`lgaussian`](https://hreinwald.github.io/drc/reference/lgaussian.md),
[`drm`](https://hreinwald.github.io/drc/reference/drm.md)

## Author

Christian Ritz

## Examples

``` r
## Simulated bell-shaped data: fit with drm() and plot
dose <- rep(c(0, 0.5, 1, 2, 4, 8, 16, 32, 64), each = 3)
mu <- gaussian()$fct(dose, matrix(c(1.2, 1, 5, 8, 1.5), length(dose), 5, byrow = TRUE))
dat <- data.frame(dose = dose, response = mu + 0.05 * sin(seq_along(dose) * 2.3))
m1 <- drm(response ~ dose, data = dat, fct = gaussian())
plot(m1, type = "all", main = "Gaussian biphasic model")


## Effect of the baseline c (fixed at 0.5) and the peak height d (fixed at 4)
m.c <- drm(response ~ dose, data = dat, fct = gaussian(fixed = c(NA, 0.5, NA, NA, NA)))
m.d <- drm(response ~ dose, data = dat, fct = gaussian(fixed = c(NA, NA, 4, NA, NA)))
plot(m1, type = "all", main = "Fixed baseline and peak height")
plot(m.c, type = "none", add = TRUE, col = "darkgreen")
plot(m.d, type = "none", add = TRUE, col = "red")
legend("topright", c("free", "c = 0.5", "d = 4"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")

```
