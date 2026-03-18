# The log-logistic function

A very general way of specifying log-logistic models under various
constraints on parameters.

## Usage

``` r
llogistic(
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

  numeric vector of length 5, specifying fixed parameters (use NA for
  non-fixed parameters).

- names:

  character vector of length 5, specifying the names of the parameters:
  b, c, d, e, f.

- method:

  character string indicating the self starter function to use.

- ssfct:

  a self starter function to be used.

- fctName:

  optional character string used internally.

- fctText:

  optional character string used internally.

## Value

A list containing the nonlinear function, the self starter function, and
the parameter names.

## Details

The five-parameter log-logistic function is given by the expression
\$\$f(x) = c + \frac{d-c}{(1+\exp(b(\log(x)-\log(e))))^f}\$\$

## References

Finney, D. J. (1979).

Seber, G. A. F. and Wild, C. J. (1989).

## See also

[`LL.2`](https://hreinwald.github.io/drc/reference/LL.2.md),
[`LL.3`](https://hreinwald.github.io/drc/reference/LL.3.md),
[`LL.4`](https://hreinwald.github.io/drc/reference/LL.4.md),
[`LL.5`](https://hreinwald.github.io/drc/reference/LL.5.md)

## Author

Christian Ritz
