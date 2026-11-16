# U-shaped CRS five-parameter model (alpha=1)

Five-parameter u-shaped CRS hormesis model with alpha=1.

## Usage

``` r
UCRS.5a(names = c("b", "c", "d", "e", "f"), ...)
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
[`CRS.5a`](https://hreinwald.github.io/drc/reference/CRS.5a.md)

## Examples

``` r
## Simulated U-shaped data (alpha = 1): fit with drm() and plot
dose <- rep(c(0.1, 0.3, 1, 3, 10, 30, 100, 300), each = 3)
mu <- ucedergreen(alpha = 1)$fct(dose, matrix(c(1.5, 1, 6, 10, -2), length(dose), 5, byrow = TRUE))
dat <- data.frame(dose = dose, response = mu + 0.15 * sin(seq_along(dose) * 2.3))
m1 <- drm(response ~ dose, data = dat, fct = UCRS.5a())
plot(m1, type = "all", main = "UCRS.5a")

```
