# Cook's distance for nonlinear dose-response models

Cook's distance values are provided for nonlinear dose-response model
fits using the same formulas as in linear regression but based on the
corresponding approximate quantities available for nonlinear models.

## Usage

``` r
# S3 method for class 'drc'
cooks.distance(model, ...)
```

## Arguments

- model:

  an object of class 'drc'.

- ...:

  additional arguments (not used).

## Value

A vector of Cook's distance values, one value per observation.

## Author

Christian Ritz

## Examples

``` r
ryegrass.LL.4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

cooks.distance(ryegrass.LL.4)
#>            1            2            3            4            5            6 
#> 7.453159e-03 7.044772e-03 4.714696e-02 4.844894e-02 2.870894e-02 4.723940e-03 
#>            7            8            9           10           11           12 
#> 6.453374e-02 4.817127e-02 3.034449e-03 1.086166e-01 1.026316e-03 1.159960e-01 
#>           13           14           15           16           17           18 
#> 6.500257e-01 1.505664e-02 6.990776e-01 8.318727e-03 1.370597e-03 1.649069e-03 
#>           19           20           21           22           23           24 
#> 3.231490e-03 6.070437e-05 1.244105e-02 1.159916e-02 1.468742e-02 4.949825e-04 
```
