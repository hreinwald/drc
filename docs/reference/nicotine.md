# nicotine

Data from an acute toxicity test with nicotine. For each of several
concentrations, the total number of subjects and the number of dead
subjects were recorded.

## Usage

``` r
data(nicotine)
```

## Format

A data frame with 12 observations on the following 3 variables.

- `conc`:

  a numeric vector

- `total`:

  a numeric vector

- `num.dead`:

  a numeric vector

## Examples

``` r
library(drc)

## Displaying the data
head(nicotine)
#>     conc total num.dead
#> 1 0.0000    45        3
#> 2 0.0025    50        5
#> 3 0.0050    46        4
#> 4 0.0100    50        3
#> 5 0.0200    46       11
#> 6 0.0300    46       20

## Fitting a two-parameter log-logistic model for binomial response
nicotine.m1 <- drm(num.dead/total ~ conc, weights = total,
data = nicotine, fct = LL.2(), type = "binomial")
summary(nicotine.m1)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) with lower limit at 0 and upper limit at 1 (2 parms)
#> 
#> Parameter estimates:
#> 
#>                 Estimate Std. Error t-value   p-value    
#> b:(Intercept) -1.7590301  0.1470173 -11.965 < 2.2e-16 ***
#> e:(Intercept)  0.0288876  0.0021099  13.691 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Plotting the fitted curve
plot(nicotine.m1, xlab = "Concentration", ylab = "Proportion dead", ylim = c(0, 1))
```
