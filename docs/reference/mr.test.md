# Mizon-Richard test for dose-response models

The function provides a lack-of-fit test for the mean structure based on
the Mizon-Richard test as compared to a specific alternative model.

## Usage

``` r
mr.test(object1, object2, object, x, var.equal = TRUE, component = 1)
```

## Arguments

- object1:

  object of class 'drc' (null model).

- object2:

  object of class 'drc' (alternative model).

- object:

  object of class 'drc' (fitted model under alternative).

- x:

  numeric vector of dose values.

- var.equal:

  logical indicating whether or not equal variances can be assumed
  across doses.

- component:

  numeric vector specifying the component(s) in the parameter vector to
  use in the test.

## Value

A p-value for test of the null hypothesis that the chosen mean structure
is appropriate as compared to the alternative mean structure provided
(see Ritz and Martinussen (2011) for a detailed explanation).

## Details

The function provides a p-value indicating whether or not the mean
structure is appropriate.

The test is applicable even in cases where data are non-normal or
exhibit variance heterogeneity.

## Note

This functionality is still experimental: Currently, the null and
alternative models are hardcoded! In the future the function will be
working for null and alternative models specified by the user.

## References

Ritz, C and Martinussen, T. (2011) Lack-of-fit tests for assessing mean
structures for continuous dose-response data, *Environmental and
Ecological Statistics*, **18**, 349–366

## See also

See also
[`modelFit`](https://hreinwald.github.io/drc/reference/modelFit.md) for
details on the related lack-of-fit test against an ANOVA model.

## Author

Christian Ritz

## Examples

``` r
## Fitting log-logistic and Weibull models
## The Weibull model is the alternative
etmotc.m1<-drm(rgr1~dose1, data=etmotc[1:15,], fct=LL.4())
etmotc.m2 <- update(etmotc.m1, fct=W1.4())

## Fitting the fitted model (using the alternative model)
etmotc.m3 <- drm(fitted(etmotc.m1)~dose1, data=etmotc[1:15,], fct=W1.4())

## Handling missing values
xVec <- etmotc[1:15,]$dose1
xVec[1:8] <- 1e-10  # avoiding 0's

## Obtaining the Mizon-Richard test
mr.test(etmotc.m1, etmotc.m2, etmotc.m3, xVec, var.equal = FALSE)
#>   Statistic     p-value  Difference          SE 
#> -1.65084985  0.09876924 -0.01936982  0.01173324 
```
