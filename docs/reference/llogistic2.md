# Five-Parameter Log-Logistic Model with log(ED50) as Parameter

A five-parameter log-logistic model where the ED50 is parameterised on
the log scale. The mean function is: \$\$f(x) = c + \frac{d - c}{(1 +
\exp(b(\log(x) - e)))^f}\$\$ where `e` is the logarithm of the ED50 (not
exponentiated).

## Usage

``` r
llogistic2(
  fixed = c(NA, NA, NA, NA, NA),
  names = c("b", "c", "d", "e", "f"),
  ss = c("1", "2", "3"),
  ssfct = NULL,
  fctName,
  fctText
)
```

## Arguments

- fixed:

  numeric vector of length 5. Specifies which parameters are fixed and
  at what value. Use `NA` for parameters that should be estimated.

- names:

  character vector of length 5 giving the names of the parameters `b`,
  `c`, `d`, `e`, and `f`.

- ss:

  character string indicating the self-starter version to use. One of
  `"1"` (default), `"2"`, or `"3"`.

- ssfct:

  optional self-starter function. If provided, overrides the built-in
  self-starter selected by `ss`.

- fctName:

  optional character string specifying the name of the function.

- fctText:

  optional character string providing a short description of the
  function.

## Value

A list of class `"llogistic"` containing the nonlinear function,
self-starter function, parameter names, and related helper functions.

## See also

[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md),
[`LL2.2`](https://hreinwald.github.io/drc/reference/LL2.2.md),
[`LL2.3`](https://hreinwald.github.io/drc/reference/LL2.3.md),
[`LL2.4`](https://hreinwald.github.io/drc/reference/LL2.4.md),
[`LL2.5`](https://hreinwald.github.io/drc/reference/LL2.5.md)

## Author

Christian Ritz
