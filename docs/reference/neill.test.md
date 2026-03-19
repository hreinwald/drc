# Neill's lack-of-fit test for dose-response models

`neill.test` provides a lack-of-fit test for non-linear regression
models. It is applicable both in cases where there are replicates (in
which case it reduces to the standard lack-of-fit test against an ANOVA
model) and in cases where there are no replicates, though then a
grouping has to be provided.

## Usage

``` r
neill.test(
  object,
  grouping,
  method = c("c-finest", "finest", "percentiles"),
  breakp = NULL,
  display = TRUE
)
```

## Arguments

- object:

  object of class 'drc' or 'nls'.

- grouping:

  character or numeric vector that provides the grouping of the dose
  values.

- method:

  character string specifying the method to be used to generate a
  grouping of the dose values.

- breakp:

  numeric vector of break points for generating dose intervals that form
  a grouping.

- display:

  logical. If TRUE results are displayed. Otherwise they are not (useful
  in simulations).

## Value

The function returns an object of class anova which is displayed using
`print.anova`.

## Details

The functions use the methods
[`df.residual`](https://rdrr.io/r/stats/df.residual.html) and
[`residuals`](https://rdrr.io/r/stats/residuals.html) and the `data`
component of `object` (only for determining the number of observations).

## Note

A clustering technique could be employed to determine the grouping to be
used in cases where there are no replicates. There should at most be
ceiling(n/2) clusters as otherwise some observations will not be used in
the test. At the other end there need to be more clusters than
parameters in the model.

## References

Neill, J. W. (1988) Testing for lack of fit in nonlinear regression,
*Ann. Statist.*, **16**, 733–740

## See also

See also
[`modelFit`](https://hreinwald.github.io/drc/reference/modelFit.md) for
details on the lack-of-fit test against an ANOVA model.

## Author

Christian Ritz

## Examples

``` r
### Example with 'drc' object

## Lack-of-fit test against ANOVA
ryegrass.m1 <-drm(rootl~conc, data = ryegrass, fct = LL.4())
modelFit(ryegrass.m1)
#> Lack-of-fit test
#> 
#>           ModelDf    RSS Df F value p value
#> ANOVA          17 5.1799                   
#> DRC model      20 5.4002  3  0.2411  0.8665

## The same test using 'neill.test'
neill.test(ryegrass.m1, ryegrass$conc)
#> Grouping used
#> 
#>    0 0.94 1.88 3.75  7.5   15   30 
#>    6    3    3    3    3    3    3 
#> 
#> Neill's lack-of-fit test
#> 
#>  F value p value
#>   0.2411  0.8665

## Generating a grouping
neill.test(ryegrass.m1, method="c-finest")
#> Grouping used
#> 
#> 1 2 3 4 5 6 7 
#> 6 3 3 3 3 3 3 
#> 
#> Neill's lack-of-fit test
#> 
#>  F value p value
#>   0.2411  0.8665
neill.test(ryegrass.m1, method="finest")
#> Grouping used
#> 
#>  1  2  3  4  5  6  7  8  9 10 11 12 
#>  2  2  2  2  2  2  2  2  2  2  2  2 
#> 
#> Neill's lack-of-fit test
#> 
#>  F value p value
#>   1.0625  0.4462
neill.test(ryegrass.m1, method="perc")
#> Grouping used
#> 
#>    (-Inf,0]    (0,1.88] (1.88,3.75]   (3.75,15]   (15, Inf] 
#>           6           6           3           6           3 
#> 
#> Neill's lack-of-fit test
#> 
#>  F value p value
#>   0.7545  0.3959
```
