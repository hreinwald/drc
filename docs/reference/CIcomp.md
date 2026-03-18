# Classical combination index for effective doses

Calculates the classical combination index for effective doses in binary
mixture experiments.

## Usage

``` r
CIcomp(mixProp, modelList, EDvec)
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

## Value

A matrix with one row per ED value. Columns contain estimated
combination indices, their standard errors and 95% confidence intervals,
p-value for testing CI=1, estimated ED values for the mixture data and
assuming concentration addition (CA) with corresponding standard errors.

## References

Martin-Betancor, K. and Ritz, C. and Fernandez-Pinas, F. and Leganes, F.
and Rodea-Palomares, I. (2015) Defining an additivity framework for
mixture research in inducible whole-cell biosensors, *Scientific
Reports* **17200**.

## See also

[`CIcompX`](https://hreinwald.github.io/drc/reference/CIcompX.md),
[`plotFACI`](https://hreinwald.github.io/drc/reference/plotFACI.md),
[`mixture`](https://hreinwald.github.io/drc/reference/mixture.md)

## Author

Christian Ritz and Ismael Rodea-Palomares

## Examples

``` r
## Fitting marginal models for the 2 pure substances
acidiq.0 <- drm(rgr ~ dose, data = subset(acidiq, pct == 999 | pct == 0), fct = LL.4())
acidiq.100 <- drm(rgr ~ dose, data = subset(acidiq, pct == 999 | pct == 100), fct = LL.4())

## Fitting model for single mixture with ratio 17:83
acidiq.17 <- drm(rgr ~ dose, data = subset(acidiq, pct == 17 | pct == 0), fct = LL.4())

## Calculation of combination indices based on ED10, ED20, ED50
CIcomp(0.17, list(acidiq.17, acidiq.0, acidiq.100), c(10, 20, 50))
#>      combInd         SE     lowCI   highCI    CAdiffp     ED.CA    SE.CA
#> 10 1.7180152 0.31407144 1.1024352 2.333595 0.02224534  76.91373 11.85583
#> 20 1.3421604 0.16702874 1.0147841 1.669537 0.04050985 140.38385 14.47436
#> 50 0.9035949 0.08440138 0.7381682 1.069022 0.25336168 382.44378 32.11935
#>      ED.mix   SE.mix
#> 10 132.1390 12.98677
#> 20 188.4176 13.13050
#> 50 345.5742 14.12771
```
