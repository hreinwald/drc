# Offsprings resulting from a toxicity test

Results from a chronic reproduction toxicity test with seven different
concentrations of waste water. The response was the number of offspring
produced by the water flea *Ceriodaphnia dubia*.

## Usage

``` r
data(C.dubia)
```

## Format

A data frame with 50 observations of the following 2 variables.

- `conc`:

  a numeric vector giving waste water in percentage

- `number`:

  a numeric vector

## Source

A. J. Bailer and J. T. Oris (1997). Estimating inhibition concentrations
for different response scales using generalized linear models.
Environmental Toxicology and Chemistry, **16**:1554–1559.

## Examples

``` r
library(drc)

## Displaying the data
head(C.dubia)
#>   conc number
#> 1    0     27
#> 2    0     30
#> 3    0     29
#> 4    0     31
#> 5    0     16
#> 6    0     15

## Fitting a three-parameter log-logistic model
C.dubia.m1 <- drm(number ~ conc, data = C.dubia, fct = LL.3())
summary(C.dubia.m1)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) with lower limit at 0 (3 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept)  18.4674   145.8889  0.1266  0.899809    
#> d:(Intercept)  28.7500     1.1123 25.8474 < 2.2e-16 ***
#> e:(Intercept)  12.0207     3.7191  3.2322  0.002247 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  7.034731 (47 degrees of freedom)

## Plotting fitted curve together with the original data
plot(C.dubia.m1, xlab = "Concentration (%)", ylab = "Number of offspring")
```
