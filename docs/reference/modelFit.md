# Assessing the model fit

Checking the fit of a dose-response model by means of formal
significance tests.

## Usage

``` r
modelFit(object, test = NULL, method = c("gof", "cum"))
```

## Arguments

- object:

  object of class 'drc'.

- test:

  character string defining the test method to apply.

- method:

  character string specifying the method to be used for assessing the
  model fit.

## Value

An object of class 'anova' which will be displayed in much the same way
as an ordinary ANOVA table.

## Details

Currently two methods are available. For continuous data the classical
lack-of-fit test is applied (Bates and Watts, 1988). The test compares
the dose-response model to a more general ANOVA model using an
approximate F-test. For quantal data the crude goodness-of-fit test
based on Pearson's statistic is used.

None of these tests are very powerful. A significant test result is more
alarming than a non-significant one.

## References

Bates, D. M. and Watts, D. G. (1988) *Nonlinear Regression Analysis and
Its Applications*, New York: Wiley & Sons (pp. 103–104).

## Author

Christian Ritz

## Examples

``` r
## Comparing the four-parameter log-logistic model
##  to a one-way ANOVA model using an approximate F test
## in other words applying a lack-of-fit test
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
modelFit(ryegrass.m1)
#> Lack-of-fit test
#> 
#>           ModelDf    RSS Df F value p value
#> ANOVA          17 5.1799                   
#> DRC model      20 6.0242  3  0.9236  0.4506
```
