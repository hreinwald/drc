# Leaf length of barley

In an experiment barley was grown in a hydroponic solution with a
herbicide.

## Usage

``` r
data(leaflength)
```

## Format

A data frame with 42 observations on the following 2 variables.

- `Dose`:

  a numeric vector

- `DW`:

  a numeric vector

## Details

The dataset exhibits a large hormetical effect.

## Source

Nina Cedergreen, Royal Veterinary and Agricultural University, Denmark.

## Examples

``` r
library(drc)

## Fitting a hormesis model
leaflength.crs4c1 <- drm(DW ~ Dose, data = leaflength, fct = CRS.4c())
plot(fitted(leaflength.crs4c1), residuals(leaflength.crs4c1))


leaflength.crs4c2 <- boxcox(drm(DW ~ Dose, data = leaflength, fct = CRS.4c()), 
method = "anova", plotit = FALSE)
summary(leaflength.crs4c2)
#> 
#> Model fitted: Cedergreen-Ritz-Streibig with lower limit 0 (alpha=0.25) (4 parms)
#> 
#> Parameter estimates:
#> 
#>                 Estimate Std. Error t-value   p-value    
#> b:(Intercept)   0.489758   0.030686 15.9603 < 2.2e-16 ***
#> d:(Intercept)  10.020054   1.590491  6.3000 2.209e-07 ***
#> e:(Intercept)   0.019138   0.028983  0.6603    0.5130    
#> f:(Intercept) 381.722590 234.463424  1.6281    0.1118    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  1.188513 (38 degrees of freedom)
#> 
#> Non-normality/heterogeneity adjustment through Box-Cox transformation
#> 
#> Estimated lambda: 0.5 
#> Confidence interval for lambda: [0.420,0.785] 
#> 

## Plottinf fitted curve and original data
plot(leaflength.crs4c2, broken = TRUE, conLevel = 0.001, type = "all", legend = FALSE, 
ylab = "Produced leaf length (cm)", xlab = "Metsulfuron-methyl (mg/l)",
main = "Hormesis: leaf length of barley")
#> Warning: "conLevel" is not a graphical parameter
#> Warning: "conLevel" is not a graphical parameter
#> Warning: "conLevel" is not a graphical parameter
#> Warning: "conLevel" is not a graphical parameter
#> Warning: "conLevel" is not a graphical parameter
```
