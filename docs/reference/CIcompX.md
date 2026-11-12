# Calculation of combination index for binary mixtures

For single mixture data, combination indices for effective doses as well
as effects may be calculated. This is an extended version of
[`CIcomp`](https://hreinwald.github.io/drc/reference/CIcomp.md).

## Usage

``` r
CIcompX(mixProp, modelList, EDvec, EDonly = FALSE)
```

## Arguments

- mixProp:

  a numeric value between 0 and 1 specifying the mixture
  proportion/ratio.

- modelList:

  a list containing 3 model fits using
  [`drm`](https://hreinwald.github.io/drc/reference/drm.md): the mixture
  model fit first, followed by the 2 pure substance model fits.

- EDvec:

  a numeric vector of effect levels (percentages between 0 and 100).

- EDonly:

  logical. If TRUE, only combination indices for effective doses are
  calculated.

## Value

A list with components `Effx`, `Effy` (unless `EDonly = TRUE`), `CAx`,
`CAy` (unless `EDonly = TRUE`), and `EDvec`.

## References

Martin-Betancor, K. and Ritz, C. and Fernandez-Pinas, F. and Leganes, F.
and Rodea-Palomares, I. (2015) Defining an additivity framework for
mixture research in inducible whole-cell biosensors, *Scientific
Reports* **17200**.

## See also

[`CIcomp`](https://hreinwald.github.io/drc/reference/CIcomp.md),
[`plotFACI`](https://hreinwald.github.io/drc/reference/plotFACI.md),
[`mixture`](https://hreinwald.github.io/drc/reference/mixture.md)

## Author

Christian Ritz and Ismael Rodea-Palomares
