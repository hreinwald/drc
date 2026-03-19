# Printing key features

`print` displays brief information on an object of class 'drc'.

## Usage

``` r
# S3 method for class 'drc'
print(x, ..., digits = max(3, getOption("digits") - 3))
```

## Arguments

- x:

  an object of class 'drc'.

- ...:

  additional arguments.

- digits:

  an integer giving the number of digits of the parameter coefficients.
  Default is 3.

## Value

The object is returned invisibly.

## Author

Christian Ritz

## Examples

``` r
## Fitting a four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~conc, data = ryegrass, fct = LL.4())

## Displaying the model fit
print(ryegrass.m1)
#> 
#> A 'drc' model.
#> 
#> Call:
#> drm(formula = rootl ~ conc, data = ryegrass, fct = LL.4())
#> 
#> Coefficients:
#> b:(Intercept)  c:(Intercept)  d:(Intercept)  e:(Intercept)  
#>        2.9822         0.4814         7.7930         3.0580  
#> 
ryegrass.m1  # gives the same output as the previous line
#> 
#> A 'drc' model.
#> 
#> Call:
#> drm(formula = rootl ~ conc, data = ryegrass, fct = LL.4())
#> 
#> Coefficients:
#> b:(Intercept)  c:(Intercept)  d:(Intercept)  e:(Intercept)  
#>        2.9822         0.4814         7.7930         3.0580  
#> 
```
