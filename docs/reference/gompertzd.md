# Derivative of the Gompertz function

`gompertzd` provides a way of specifying the derivative of the Gompertz
function as a dose-response model.

## Usage

``` r
gompertzd(fixed = c(NA, NA), names = c("a", "b"))
```

## Arguments

- fixed:

  numeric vector. Specifies which parameters are fixed and at what value
  they are fixed. NAs for parameters that are not fixed.

- names:

  a vector of character strings giving the names of the parameters
  (should not contain ":"). The default is (notice the order): a, b.

## Value

A list containing the model function, the self starter function and the
parameter names.

## Details

The derivative of the Gompertz function is defined as \$\$f(x) = a
\exp(bx-a/b(\exp(bx)-1))\$\$ For \\a\>0\\ and \\b\\ not 0, the function
is decreasing, equaling \\a\\ at \\x=0\\ and approaching 0 at plus
infinity.

## See also

[`gompertz`](https://hreinwald.github.io/drc/reference/gompertz.md),
[`drm`](https://hreinwald.github.io/drc/reference/drm.md)

## Author

Christian Ritz
