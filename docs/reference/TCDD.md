# Liver tumor incidence

Liver tumor incidence in Sprague-Dawley rats exposed to the chemical
like 2,3,7,8-tetrachlorodibenzo-pdioxin (TCDD).

## Usage

``` r
data(TCDD)
```

## Format

A data frame with 6 observations on the following 3 variables.

- `conc`:

  a numeric vector reporting the concentration of TCDD (ng/kg)

- `total`:

  a numeric vector

- `incidence`:

  a numeric vector

## Source

R. Kociba, D. Keyes, J. Beyer, R. Carreon, C. Wade, D. Dittenber, R.
Kalnins, L. Frauson, C. Park, S. Barnard, R. Hummel, and C. Humiston
(1978). Results of a two-year chronic toxicity and oncogenicity study of
2,3,7,8-tetrachlorodibenzo-p-dioxin in rats. Toxicology and Applied
Pharmacology, **46(2)**:279–303.

## Examples

``` r
library(drc)

## Displaying the data
head(TCDD)
#>    conc total incidence
#> 1  0.00    86         2
#> 2  1.55    50         1
#> 3  7.15    50         9
#> 4 38.56    45        14

## Fitting a two-parameter log-logistic model for binomial response
TCDD.m1 <- drm(incidence/total ~ conc, weights = total,
data = TCDD, fct = LL.2(), type = "binomial")
summary(TCDD.m1)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) with lower limit at 0 and upper limit at 1 (2 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) -0.73196    0.20589 -3.5551 0.0003778 ***
#> e:(Intercept) 96.49945   58.87278  1.6391 0.1011886    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Plotting the fitted curve
plot(TCDD.m1, xlab = "Concentration of TCDD (ng/kg)", ylab = "Tumor incidence")
```
