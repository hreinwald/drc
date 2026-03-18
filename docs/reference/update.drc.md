# Updating and re-fitting a model

`update` updates and re-fits a model on the basis of an object of class
'drc'.

## Usage

``` r
# S3 method for class 'drc'
update(object, ..., evaluate = TRUE)
```

## Arguments

- object:

  an object of class 'drc'.

- ...:

  arguments to alter in object.

- evaluate:

  logical. If TRUE model is re-fit; otherwise an unevaluated call is
  returned.

## Value

An object of class 'drc'.

## Author

Christian Ritz

## Examples

``` r
## Fitting a four-parameter Weibull model
model1 <- drm(ryegrass, fct = W1.4())

## Updating 'model1' by fitting a three-parameter Weibull model instead
model2 <- update(model1, fct = W1.3())
anova(model2, model1)
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
```
