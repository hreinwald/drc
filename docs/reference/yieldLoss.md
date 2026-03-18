# Calculating yield loss parameters

Calculation of parameters in the re-parameterization of the
Michaelis-Menten model that is commonly used to assess yield loss (the
rectangular hyperbola model).

## Usage

``` r
yieldLoss(object, interval = c("none", "as"), level = 0.95, display = TRUE)
```

## Arguments

- object:

  object of class 'drc'.

- interval:

  character string specifying the type of confidence intervals. The
  default is "none". Use "as" for asymptotically-based confidence
  intervals.

- level:

  numeric. The level for the confidence intervals. The default is 0.95.

- display:

  logical. If TRUE results are displayed. Otherwise they are not (useful
  in simulations).

## Value

For each of the two parameters, a matrix with two or more columns,
containing the estimates and the corresponding estimated standard errors
and possibly lower and upper confidence limits.

## Details

The rectangular hyperbola model is a reparameterization of the
Michaelis-Menten in terms of parameters \\A\\ and \\I\\: \$\$Y_L =
\frac{Id}{1+Id/A}\$\$ where \\d\\ denotes the weed density and \\Y_L\\
the resulting yield loss.

## Note

This function is only for use with model fits based on Michaelis-Menten
models.

## References

Cousens, R. (1985). A simple model relating yield loss to weed density,
*Ann. Appl. Biol.*, **107**, 239–252.

## Author

Christian Ritz

## Examples

``` r
## Fitting Michaelis-Menten model
met.mm.m1 <- drm(gain~dose, product, data = methionine, fct = MM.3(),
pmodels = list(~1, ~factor(product), ~factor(product)))
#> Control measurements detected for level: control

## Yield loss parameters with standard errors
yieldLoss(met.mm.m1)
#> 
#> Estimated A parameters
#> 
#>     Estimate Std. Error
#> DLM 1736.141     18.922
#> MHA 1868.517     43.930
#> 
#> 
#> Estimated I parameters
#> 
#>     Estimate Std. Error
#> DLM  44578.0    11225.6
#> MHA  16827.3     3942.7

## Also showing confidence intervals
yieldLoss(met.mm.m1, "as")
#> 
#> Estimated A parameters
#> 
#>     Estimate Std. Error    Lower    Upper
#> DLM 1736.141     18.922 1683.606 1788.676
#> MHA 1868.517     43.930 1746.547 1990.487
#> 
#> 
#> Estimated I parameters
#> 
#>     Estimate Std. Error   Lower   Upper
#> DLM  44578.0    11225.6 13410.7 75745.2
#> MHA  16827.3     3942.7  5880.7 27773.8
```
