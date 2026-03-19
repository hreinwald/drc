# Five-parameter Weibull type 2 model with lag time

A five-parameter Weibull type 2 model extended with a lag time parameter
`t0`. The model is given by the expression \$\$f(x) = c + (d - c)(1 -
\exp(-\exp(b(\log(x - t0) - \log(e)))))\$\$ for \\x \> t0\\ and \\f(x) =
c\\ otherwise.

## Usage

``` r
weibull2x(
  fixed = rep(NA, 5),
  names = c("b", "c", "d", "e", "t0"),
  method = c("1", "2", "3", "4"),
  ssfct = NULL,
  fctName,
  fctText
)
```

## Arguments

- fixed:

  numeric vector of length 5. Specifies which parameters are fixed and
  at what value. Use `NA` for parameters that should be estimated
  (default is `rep(NA, 5)`).

- names:

  character vector of length 5 giving the names of the parameters
  (default is `c("b", "c", "d", "e", "t0")`).

- method:

  character string indicating the self starter method to use. One of
  `"1"`, `"2"`, `"3"`, or `"4"`.

- ssfct:

  a self starter function. If `NULL` (default), a built-in self starter
  is used.

- fctName:

  optional character string specifying the function name (used
  internally).

- fctText:

  optional character string specifying the function description (used
  internally).

## Value

A list of class `"Weibull-2"` containing the nonlinear function, self
starter function, and parameter names.

## Details

The lag time parameter `t0` cannot be fixed.

## See also

[`weibull2`](https://hreinwald.github.io/drc/reference/weibull2.md),
[`W2x.3`](https://hreinwald.github.io/drc/reference/W2x.3.md),
[`W2x.4`](https://hreinwald.github.io/drc/reference/W2x.4.md)

## Author

Christian Ritz
