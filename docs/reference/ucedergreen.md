# U-shaped Cedergreen-Ritz-Streibig model

`ucedergreen` provides a very general way of specifying the
Cedergreen-Ritz-Streibig modified log-logistic model for describing
u-shaped hormesis, under various constraints on the parameters.

## Usage

``` r
ucedergreen(
  fixed = c(NA, NA, NA, NA, NA),
  names = c("b", "c", "d", "e", "f"),
  method = c("1", "2", "3", "4"),
  ssfct = NULL,
  alpha
)
```

## Arguments

- fixed:

  numeric vector. Specifies which parameters are fixed and at what value
  they are fixed. NAs for parameters that are not fixed.

- names:

  a vector of character strings giving the names of the parameters
  (should not contain ":"). The order of the parameters is: b, c, d, e,
  f.

- method:

  character string indicating the self starter function to use.

- ssfct:

  a self starter function to be used.

- alpha:

  numeric value between 0 and 1, reflecting the steepness of the
  hormesis peak. This argument must be specified.

## Value

A list containing the non-linear function, the self starter function and
the parameter names.

## Details

The u-shaped model is given by the expression \$\$f(x) = c + d -
\frac{d-c+f \exp(-1/x^{\alpha})}{1+\exp(b(\log(x)-\log(e)))}\$\$

## References

Cedergreen, N. and Ritz, C. and Streibig, J. C. (2005) Improved
empirical models describing hormesis, *Environmental Toxicology and
Chemistry* **24**, 3166–3172.

## See also

[`cedergreen`](https://hreinwald.github.io/drc/reference/cedergreen.md),
[`UCRS.4a`](https://hreinwald.github.io/drc/reference/UCRS.4a.md),
[`UCRS.5a`](https://hreinwald.github.io/drc/reference/UCRS.5a.md),
[`drm`](https://hreinwald.github.io/drc/reference/drm.md)

## Author

Christian Ritz
