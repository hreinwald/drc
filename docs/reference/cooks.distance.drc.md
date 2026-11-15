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
#> 0.0074533370 0.0070445992 0.0471465147 0.0484493866 0.0287092899 0.0047237984 
#>            7            8            9           10           11           12 
#> 0.0645333094 0.0481716296 0.0030343572 0.1086155004 0.0010262102 0.1159971709 
#>           13           14           15           16           17           18 
#> 0.6500270558 0.0150568459 0.6990762170 0.0083187135 0.0013706026 0.0016490630 
#>           19           20           21           22           23           24 
#> 0.0032314889 0.0000607045 0.0124410494 0.0115991555 0.0146874200 0.0004949825 
```
