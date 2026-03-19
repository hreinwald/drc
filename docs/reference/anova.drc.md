# ANOVA for dose-response model fits

`anova` produces an analysis of variance table for one or two non-linear
model fits.

## Usage

``` r
# S3 method for class 'drc'
anova(object, ..., details = TRUE, test = NULL)
```

## Arguments

- object:

  an object of class 'drc'.

- ...:

  additional arguments.

- details:

  logical indicating whether or not details on the models compared
  should be displayed. Default is TRUE (details are displayed).

- test:

  a character string specifying the test statistic to be applied. Use
  `"od"` to assess overdispersion for binomial data.

## Value

An object of class 'anova'.

## Details

Specifying only a single object gives a test for lack-of-fit, comparing
the non-linear regression model to a more general one-way or two-way
ANOVA model.

If two objects are specified a test for reduction from the larger to the
smaller model is given. (This only makes statistical sense if the models
are nested, that is: one model is a submodel of the other model.)

## Author

Christian Ritz

## Examples

``` r
## Comparing Weibull three- and four-parameter models using an F test
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = W1.4())
ryegrass.m2 <- drm(rootl ~ conc, data = ryegrass, fct = W1.3())
anova(ryegrass.m2, ryegrass.m1)
#> 
#> 1st model
#>  fct:      W1.3()
#> 2nd model
#>  fct:      W1.4()
#> 
#> ANOVA table
#> 
#>           ModelDf    RSS Df F value p value
#> 1st model      21 8.9520                   
#> 2nd model      20 6.0242  1  9.7205  0.0054

anova(ryegrass.m2, ryegrass.m1, details = FALSE)  # without details
#> ANOVA table
#> 
#>           ModelDf    RSS Df F value p value
#> 1st model      21 8.9520                   
#> 2nd model      20 6.0242  1  9.7205  0.0054
```
