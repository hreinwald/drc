# Calculation of backfit values from a fitted dose-response model

By inverse regression backfitted dose values are calculated for the mean
response per dose.

## Usage

``` r
backfit(drcObject)
```

## Arguments

- drcObject:

  an object of class 'drc'.

## Value

Two columns with the original dose values and the corresponding
backfitted values using the fitted dose-response model. For extreme dose
values (e.g., high dose) the backfitted values may not be well-defined.

## See also

A related function is
[`ED.drc`](https://hreinwald.github.io/drc/reference/ED.drc.md).

## Author

Christian Ritz after a suggestion from Keld Sorensen.

## Examples

``` r
ryegrass.LL.4 <- drm(rootl~conc, data=ryegrass, fct=LL.4())

backfit(ryegrass.LL.4)
#> Warning: NaNs produced
#>       dose   Estimate
#> [1,]  0.00  0.5500692
#> [2,]  0.94  0.7743783
#> [3,]  1.88  1.8744292
#> [4,]  3.75  3.7830500
#> [5,]  7.50  7.0811832
#> [6,] 15.00 10.1667582
#> [7,] 30.00        NaN
```
