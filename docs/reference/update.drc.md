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

## Details

The arguments stored in the call of the fit (e.g. `formula = myFormula`,
`data = myData`) are evaluated in the environment in which
[`drm()`](https://hreinwald.github.io/drc/reference/drm.md) was called,
so [`update()`](https://rdrr.io/r/stats/update.html) - and functions
that refit models through it, such as
[`mselect`](https://hreinwald.github.io/drc/reference/mselect.md) and
[`maED`](https://hreinwald.github.io/drc/reference/maED.md) - also work
for fits created inside functions or in a separate environment, and
never pick up an object of the same name from the global environment by
accident. New arguments supplied in `...` are evaluated in the
environment [`update()`](https://rdrr.io/r/stats/update.html) is called
from.

If `data` is not supplied and the original `data` expression can no
longer be found or now refers to different data, the data stored in the
fit (`object$origData`) is used (with a warning in the latter case).

Fits made inside a function keep a reference to that function's
environment (as [`lm()`](https://rdrr.io/r/stats/lm.html) does through
its formula), which is also saved by
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html).

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
