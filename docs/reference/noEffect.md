# Testing if there is a dose effect at all

A significance test is provided for the comparison of the dose-response
model considered and the simple linear regression model with slope 0 (a
horizontal regression line corresponding to no dose effect).

## Usage

``` r
noEffect(object)
```

## Arguments

- object:

  an object of class 'drc'.

## Value

The likelihood ratio test statistic and the corresponding degrees of
freedom and p-value are reported.

## Details

Perhaps useful for screening purposes. The test is a likelihood ratio
test against the intercept-only model and is available for the response
types `"continuous"`, `"binomial"` and `"Poisson"`; for other types
(`"negbin1"`, `"negbin2"`, `"event"`, `"ssd"`) an error is raised. A
warning is given if the dose-response model has no more free parameters
than the null model (e.g. when most parameters are fixed).

## Author

Christian Ritz

## Examples

``` r
ryegrass.LL.4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

noEffect(ryegrass.LL.4)
#> Chi-square test              Df         p-value 
#>    9.187776e+01    3.000000e+00    8.653580e-20 
# p-value < 0.0001: there is a highly significant dose effect!
```
