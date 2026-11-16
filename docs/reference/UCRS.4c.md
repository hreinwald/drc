# U-shaped CRS model with lower limit 0 (alpha=0.25)

Four-parameter u-shaped CRS hormesis model with lower limit fixed at 0
and alpha=0.25.

## Usage

``` r
UCRS.4c(names = c("b", "d", "e", "f"), ...)
```

## Arguments

- names:

  a vector of character strings giving the names of the parameters.

- ...:

  additional arguments passed to
  [`ucedergreen`](https://hreinwald.github.io/drc/reference/ucedergreen.md).

## Value

A list (see
[`ucedergreen`](https://hreinwald.github.io/drc/reference/ucedergreen.md)).

## See also

[`ucedergreen`](https://hreinwald.github.io/drc/reference/ucedergreen.md),
[`UCRS.4a`](https://hreinwald.github.io/drc/reference/UCRS.4a.md),
[`CRS.4c`](https://hreinwald.github.io/drc/reference/CRS.4c.md)

## Examples

``` r
## Simulated U-shaped data (c = 0, alpha = 0.25): fit with drm() and plot
dose <- rep(c(0.1, 0.3, 1, 3, 10, 30, 100, 300), each = 3)
mu <- ucedergreen(alpha = 0.25, fixed = c(NA, 0, NA, NA, NA))$fct(dose,
  matrix(c(1.5, 6, 10, -2), length(dose), 4, byrow = TRUE))
dat <- data.frame(dose = dose, response = mu + 0.15 * sin(seq_along(dose) * 2.3))
m1 <- drm(response ~ dose, data = dat, fct = UCRS.4c())
plot(m1, type = "all", main = "UCRS.4c")

```
