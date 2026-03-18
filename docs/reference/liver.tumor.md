# Liver tumor incidence

Liver tumor incidence in female Sprague-Dawley rats exposed to the
chemical like 2,3,7,8-tetrachlorodibenzo-pdioxin (TCDD).

## Usage

``` r
data(liver.tumor)
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

National Toxicology Program. NTP technical report on the toxicology and
carcinogenesis studies of 2,3,7,8-tetrachlorodibenzo-p-dioxin (tcdd)
(CAS No. 1746-01-6) in female Harlan Sprague-Dawley rats (gavage
studies). National Toxicology Program technical report series,
(521):4–232, apr 2006.

## Examples

``` r
library(drc)

## Displaying the data
head(liver.tumor)
#>    conc total incidence
#> 1  0.00    49         0
#> 2  2.56    48         0
#> 3  5.69    46         0
#> 4  9.79    50         0
#> 5 16.57    49         1
#> 6 29.70    53        13

## Fitting a two-parameter log-logistic model for binomial response
liver.tumor.m1 <- drm(incidence/total ~ conc, weights = total,
data = liver.tumor, fct = LL.2(), type = "binomial")
summary(liver.tumor.m1)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) with lower limit at 0 and upper limit at 1 (2 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept)  -4.9750     1.6897 -2.9443  0.003237 ** 
#> e:(Intercept)  37.1774     4.1838  8.8859 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Plotting the fitted curve
plot(liver.tumor.m1, xlab = "Concentration of TCDD (ng/kg)",
ylab = "Tumor incidence")
```
