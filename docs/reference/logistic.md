# The general asymmetric five-parameter logistic model

The five-parameter logistic model given by the expression \$\$f(x) = c +
\frac{d - c}{(1 + \exp(b(x - e)))^f}\$\$

## Usage

``` r
logistic(
  fixed = c(NA, NA, NA, NA, NA),
  names = c("b", "c", "d", "e", "f"),
  method = c("1", "2", "3", "4"),
  ssfct = NULL,
  fctName,
  fctText
)
```

## Arguments

- fixed:

  numeric vector of length 5. Specifies which parameters are fixed and
  at what value they are fixed. `NA` indicates that the corresponding
  parameter is not fixed.

- names:

  character vector of length 5 giving the names of the parameters
  `(b, c, d, e, f)`. Default is `c("b", "c", "d", "e", "f")`.

- method:

  character string indicating the self starter function to use (`"1"`,
  `"2"`, `"3"`, or `"4"`).

- ssfct:

  a self starter function to be used. If `NULL` (default), a built-in
  self starter is selected via `method`.

- fctName:

  optional character string used internally to overwrite the function
  name.

- fctText:

  optional character string used internally to overwrite the description
  text.

## Value

A list of class `"Boltzmann"` containing the nonlinear function, self
starter function, and parameter names.

## Details

This model differs from the log-logistic in that it uses `x` directly
rather than `log(x)`. It is sometimes referred to as the Boltzmann
model.

## See also

[`L.3`](https://hreinwald.github.io/drc/reference/L.3.md),
[`L.4`](https://hreinwald.github.io/drc/reference/L.4.md),
[`L.5`](https://hreinwald.github.io/drc/reference/L.5.md),
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md)

## Author

Christian Ritz

## Examples

``` r
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = L.4())
```
