# Ryegrass

Data from a dose-response experiment with ryegrass (*Lolium* sp.).
Biomass was measured at different dose levels and at two time points.

## Usage

``` r
data(ryegrass2)
```

## Format

A data frame with 27 observations on the following 3 variables.

- `dose`:

  a numeric vector

- `biomass`:

  a numeric vector

- `day`:

  a categorial vector

## Examples

``` r
library(drc)

## Displaying the data
head(ryegrass2)
#>   dose biomass day
#> 1    0    77.5   0
#> 2    0    78.0   0
#> 3    0    75.0   0
#> 4    0   214.4  15
#> 5    0   215.4  15
#> 6    0   227.5  15

## Fitting a four-parameter log-logistic model with separate curves per day
ryegrass2.m1 <- drm(biomass ~ dose, day, data = ryegrass2, fct = LL.4())
#> Control measurements detected for level: 0
summary(ryegrass2.m1)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) (4 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept)   2.9412     2.4262  1.2123  0.237709    
#> c:(Intercept)  82.7372    21.4655  3.8544  0.000807 ***
#> d:(Intercept) 179.8062    12.4500 14.4423 5.046e-13 ***
#> e:(Intercept)  13.6913     6.2483  2.1912  0.038827 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  44.3814 (23 degrees of freedom)

## Plotting the fitted curves
plot(ryegrass2.m1, xlab = "Dose", ylab = "Biomass")
```
