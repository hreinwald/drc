# Two-parameter exponential decay model

A two-parameter exponential decay model with the slope parameter `b`
fixed at 1 and the lower limit fixed at 0.

## Usage

``` r
EXD.2(fixed = c(NA, NA), names = c("d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 2. Specifies which parameters are fixed and
  at what value. Use `NA` for parameters that are not fixed.

- names:

  character vector of length 2 giving the names of the parameters. The
  default is `c("d", "e")`.

- ...:

  additional arguments passed to
  [`weibull1`](https://hreinwald.github.io/drc/reference/weibull1.md),
  most notably `method` (a character string: `"1"` (default), `"2"`,
  `"3"`, or `"4"`) which selects the self-starter method for obtaining
  starting values. See
  [`weibull1`](https://hreinwald.github.io/drc/reference/weibull1.md)
  for details.

## Value

A list of class `Weibull-1` containing the nonlinear function, self
starter function, and parameter names.

## Details

The model is given by the expression \$\$f(x) = d \exp(-x/e)\$\$

This is a special case of the Weibull type 1 model
([`weibull1`](https://hreinwald.github.io/drc/reference/weibull1.md))
with the slope fixed at 1 and the lower limit fixed at 0.

## References

Seber, G. A. F. and Wild, C. J. (1989) *Nonlinear Regression*, New York:
Wiley & Sons (pp. 338–339).

## See also

[`EXD.3`](https://hreinwald.github.io/drc/reference/EXD.3.md),
[`AR.2`](https://hreinwald.github.io/drc/reference/AR.2.md),
[`AR.3`](https://hreinwald.github.io/drc/reference/AR.3.md),
[`weibull1`](https://hreinwald.github.io/drc/reference/weibull1.md)

## Examples

``` r
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = EXD.2())
```
