# Five-parameter generalized logistic model

A five-parameter generalized logistic model (asymmetric when `f != 1`),
given by \$\$f(x) = c + \frac{d - c}{(1 + \exp(b(x - e)))^f}\$\$

## Usage

``` r
L.5(fixed = c(NA, NA, NA, NA, NA), names = c("b", "c", "d", "e", "f"), ...)
```

## Arguments

- fixed:

  numeric vector of length 5. Specifies which parameters are fixed and
  at what value they are fixed. `NA` indicates that the corresponding
  parameter is not fixed.

- names:

  character vector of length 5 giving the names of the parameters
  `(b, c, d, e, f)`. Default is `c("b", "c", "d", "e", "f")`.

- ...:

  additional arguments passed to
  [`logistic`](https://hreinwald.github.io/drc/reference/logistic.md).

## Value

A list of class `"Boltzmann"` containing the nonlinear function, self
starter function, and parameter names.

## See also

[`logistic`](https://hreinwald.github.io/drc/reference/logistic.md),
[`L.3`](https://hreinwald.github.io/drc/reference/L.3.md),
[`L.4`](https://hreinwald.github.io/drc/reference/L.4.md)

## Examples

``` r
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = L.5())
```
