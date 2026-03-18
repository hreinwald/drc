# Three-parameter logistic model

A three-parameter logistic model with the lower limit fixed at 0, given
by \$\$f(x) = \frac{d}{1 + \exp(b(x - e))}\$\$

## Usage

``` r
L.3(fixed = c(NA, NA, NA), names = c("b", "d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 3. Specifies which parameters are fixed and
  at what value they are fixed. `NA` indicates that the corresponding
  parameter is not fixed.

- names:

  character vector of length 3 giving the names of the parameters
  `(b, d, e)`. Default is `c("b", "d", "e")`.

- ...:

  additional arguments passed to
  [`logistic`](https://hreinwald.github.io/drc/reference/logistic.md).

## Value

A list of class `"Boltzmann"` containing the nonlinear function, self
starter function, and parameter names.

## See also

[`logistic`](https://hreinwald.github.io/drc/reference/logistic.md),
[`L.4`](https://hreinwald.github.io/drc/reference/L.4.md),
[`L.5`](https://hreinwald.github.io/drc/reference/L.5.md)

## Examples

``` r
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = L.3())
```
