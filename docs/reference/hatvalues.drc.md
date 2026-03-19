# Model diagnostics for nonlinear dose-response models

Hat values (leverage values) are provided for nonlinear dose-response
model fits using the same formulas as in linear regression but based on
the corresponding approximate quantities available for nonlinear models.

## Usage

``` r
# S3 method for class 'drc'
hatvalues(model, ...)
```

## Arguments

- model:

  an object of class 'drc'.

- ...:

  additional arguments (not used).

## Value

A vector of leverage values (hat values), one value per observation.

## Details

Hat values are calculated using the formula given by Cook et al. (1986)
and McCullagh and Nelder (1989). The output values can be assessed in
the same way as in linear regression.

## References

Cook, R. D. and Tsai, C.-L. and Wei, B. C. (1986) Bias in Nonlinear
Regression, *Biometrika* **73**, 615–623.

McCullagh, P. and Nelder, J. A. (1989) *Generalized Linear Models*,
Second edition, Chapman & Hall/CRC.

## Author

Christian Ritz

## Examples

``` r
ryegrass.LL.4 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

hatvalues(ryegrass.LL.4)
#>          1          2          3          4          5          6          7 
#> 0.13332342 0.13332342 0.13332342 0.13332342 0.13332342 0.13332342 0.09539668 
#>          8          9         10         11         12         13         14 
#> 0.09539668 0.09539668 0.27948291 0.27948291 0.27948291 0.28191831 0.28191831 
#>         15         16         17         18         19         20         21 
#> 0.28191831 0.12467564 0.12467564 0.12467564 0.12949336 0.12949336 0.12949336 
#>         22         23         24 
#> 0.15571960 0.15571960 0.15571960 
```
