# Measurements of germination tubes for Giant Kelp

Giant kelp, *Macrocystis pyrifera*, was exposed to 8 different
concentrations of copper and the response measured was the length of the
germination tube.

## Usage

``` r
data(GiantKelp)
```

## Format

A data frame with 39 observations of the following 2 variables.

- `dose`:

  a numeric vector

- `tubeLength`:

  a numeric vector giving the length of the germination tube (mm)

## Source

G. A. Chapman, D. L. Denton, and J. M. Lazorchak (1995). Short-term
methods for estimating the chronic toxicity of effluents and receiving
waters to west coast marine and estuarine organisms.

## Examples

``` r
library(drc)

## Displaying the data
head(GiantKelp)
#>   tubeLength dose
#> 1      19.58  0.0
#> 2      18.75  0.0
#> 3      19.14  0.0
#> 4      16.50  0.0
#> 5      17.93  0.0
#> 6      18.26  5.6

## Fitting a four-parameter log-logistic model
GiantKelp.m1 <- drm(tubeLength ~ dose, data = GiantKelp, fct = LL.4())
summary(GiantKelp.m1)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) (4 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value p-value    
#> b:(Intercept)  1.19476    0.45955  2.5999 0.01357 *  
#> c:(Intercept)  4.46327    3.46849  1.2868 0.20661    
#> d:(Intercept) 18.08505    0.77922 23.2090 < 2e-16 ***
#> e:(Intercept) 53.86481   25.46202  2.1155 0.04158 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  1.688897 (35 degrees of freedom)

## Plotting the fitted curve
plot(GiantKelp.m1, xlab = "Copper concentration", ylab = "Tube length (mm)")
```
