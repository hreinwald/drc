# Effect of cadmium on growth of green alga

Green alga (*Selenastrum capricornutum*) was exposed to cadmium chloride
concentrations ranging from 5 to 80 micro g/L in geometric progression
in 4-day population growth test.

## Usage

``` r
data(S.capricornutum)
```

## Format

A data frame with 18 observations on the following 2 variables.

- `conc`:

  a numeric vector of cadmium chloride concentrations (micro g/L)

- `count`:

  a numeric vector of algal counts (10000 x cells /ml)

## Details

The data are analysed in Bruce and Versteeg (1992) using a log-normal
dose-response model (using the logarithm with base 10).

## Source

Bruce, R. D. and Versteeg, D. J. (1992) A statistical procedure for
modeling continuous toxicity data, *Environ. Toxicol. Chem.*, **11**,
1485–1494.

## Examples

``` r
library(drc)

## Fitting 3-parameter log-normal model
s.cap.m1 <- drm(count ~ conc, data = S.capricornutum, fct = LN.3())

## Residual plot
plot(fitted(s.cap.m1), residuals(s.cap.m1))


## Fitting model with transform-both-sides approach
s.cap.m2 <- boxcox(s.cap.m1, method = "anova")

summary(s.cap.m2)
#> 
#> Model fitted: Log-normal with lower limit at 0 (3 parms)
#> 
#> Parameter estimates:
#> 
#>                 Estimate Std. Error t-value   p-value    
#> b:(Intercept)  -1.000987   0.044845 -22.321 6.394e-13 ***
#> d:(Intercept) 132.078791   7.553881  17.485 2.190e-11 ***
#> e:(Intercept)  12.428206   1.100904  11.289 9.913e-09 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1551479 (15 degrees of freedom)
#> 
#> Non-normality/heterogeneity adjustment through Box-Cox transformation
#> 
#> Estimated lambda: 0.0606 
#> Confidence interval for lambda: [-0.220, 0.414] 
#> 

## Residual plot after transformation (looks better)
plot(fitted(s.cap.m2), residuals(s.cap.m2))


## Calculating ED values on log scale
ED(s.cap.m2, c(10, 20, 50), interval="delta")
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error    Lower    Upper
#> e:10  3.45451    0.49164  2.40660  4.50242
#> e:20  5.36114    0.66213  3.94985  6.77244
#> e:50 12.42821    1.10090 10.08168 14.77473

## Fitting model with ED50 as parameter
## (for comparison)
s.cap.m3 <- drm(count ~ conc, data = S.capricornutum, fct = LN.3(loge=TRUE))
s.cap.m4 <- boxcox(s.cap.m3, method = "anova")

summary(s.cap.m4)
#> 
#> Model fitted: Log-normal with lower limit at 0 (3 parms)
#> 
#> Parameter estimates:
#> 
#>                 Estimate Std. Error t-value   p-value    
#> b:(Intercept)  -1.000990   0.044846 -22.320 6.395e-13 ***
#> d:(Intercept) 132.078289   7.553924  17.485 2.191e-11 ***
#> e:(Intercept)   2.519976   0.088582  28.448 1.821e-14 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1551479 (15 degrees of freedom)
#> 
#> Non-normality/heterogeneity adjustment through Box-Cox transformation
#> 
#> Estimated lambda: 0.0606 
#> Confidence interval for lambda: [-0.220, 0.414] 
#> 
ED(s.cap.m4, c(10, 20, 50), interval = "fls")
#> 
#> Estimated effective doses
#> 
#>       Estimate Std. Error     Lower     Upper
#> e:10  3.454549   0.142321  2.550634  4.678803
#> e:20  5.361196   0.123507  4.120349  6.975726
#> e:50 12.428297   0.088582 10.289944 15.011021
```
