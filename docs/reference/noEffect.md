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

Perhaps useful for screening purposes.

## Author

Christian Ritz

## Examples

``` r
ryegrass.LL.4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

noEffect(ryegrass.LL.4)
#> Chi-square test              Df         p-value 
#>        91.87776         3.00000         0.00000 
# p-value < 0.0001: there is a highly significant dose effect!
```
