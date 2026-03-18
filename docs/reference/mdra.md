# 3T3 mouse fibroblasts and NRU assay

The toxicity of sodium valproate was tested, using the 3T3 mouse
fibroblasts and neutral red uptake (NRU) assay. 22 different experiments
were performed independently in six laboratories, using eight
concentration levels, each with six replicates on a 96-well plate. In
addition, twelve measurements were taken for the solvent control.

## Usage

``` r
data("mdra")
```

## Format

A data frame with 1320 observations on the following 4 variables.

- `LabID`:

  a factor with levels `A` `B` `C` `D` `E` `F`

- `ExperimentID`:

  a factor with levels `1` `2` `3` `4` `5` `6` `7` `8` `9` `10` `11`
  `12` `13` `14` `15` `16` `17` `18` `19` `20` `21` `22`

- `Concentration`:

  a numeric vector

- `Response`:

  a numeric vector

## Source

http://biostatistics.dkfz.de/download/mdra/MDRA_ExampleData.csv

## References

Clothier, R., Gomez-Lechon, M. J., Kinsner-Ovaskainen, A.,
Kopp-Schneider, A., O'Connor, J. E., Prieto, P., and Stanzel, S. (2013).
Comparative analysis of eight cytotoxicity assays evaluated within the
ACuteTox Project. Toxicology in vitro, 27(4):1347–1356.

## Examples

``` r
data(mdra)

## Fit a three-parameter log-logistic model
mdra.m1 <- drm(Response ~ Concentration, data = mdra, fct = LL.3())
summary(mdra.m1)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) with lower limit at 0 (3 parms)
#> 
#> Parameter estimates:
#> 
#>                 Estimate Std. Error t-value   p-value    
#> b:(Intercept) 0.87155467 0.05369552  16.231 < 2.2e-16 ***
#> d:(Intercept) 1.03125713 0.00904617 113.999 < 2.2e-16 ***
#> e:(Intercept) 0.00381995 0.00021312  17.924 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.2333374 (1317 degrees of freedom)
plot(mdra.m1, main = "MDRA dose-response")
```
