# The effect of terbuthylazin on growth rate

Test on the effect of terbuthylazin on *Lemna minor*, performed on an
aseptic culture according to the OECD-guidelines.

## Usage

``` r
data(terbuthylazin)
```

## Format

A data frame with 30 observations on the following 2 variables.

- dose:

  a numeric vector of dose values.

- rgr:

  a numeric vector of relative growth rates.

## Details

Dose is \$\$\mu l^{-1}\$\$ and rgr is the relative growth rate of
*Lemna*.

## Source

Cedergreen N. (2004). Unpublished bioassay data.

## Examples

``` r
library(drc)

## displaying first 6 rows of the data set
head(terbuthylazin)
#>   dose       rgr
#> 1    0 0.3017731
#> 2    0 0.2760291
#> 3    0 0.3145257
#> 4    0 0.2663174
#> 5    0 0.2871303
#> 6    0 0.3805772

## Fitting log-logistic model
terbuthylazin.m1 <- drm(rgr~dose, data = terbuthylazin, fct = LL.4())
summary(terbuthylazin.m1)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) (4 parms)
#> 
#> Parameter estimates:
#> 
#>                 Estimate Std. Error t-value   p-value    
#> b:(Intercept) 1.2474e+00 2.3543e-01  5.2984 1.532e-05 ***
#> c:(Intercept) 8.0293e-04 2.4913e-02  0.0322    0.9745    
#> d:(Intercept) 3.0695e-01 9.6088e-03 31.9441 < 2.2e-16 ***
#> e:(Intercept) 1.8914e+02 3.7726e+01  5.0136 3.242e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.02699266 (26 degrees of freedom)

## Fitting log-logistic model
##  with Box-Cox transformation
terbuthylazin.m2 <- boxcox(terbuthylazin.m1, method = "anova")

summary(terbuthylazin.m2)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) (4 parms)
#> 
#> Parameter estimates:
#> 
#>                 Estimate Std. Error t-value   p-value    
#> b:(Intercept) 1.3226e+00 2.3707e-01  5.5788 7.346e-06 ***
#> c:(Intercept) 5.6549e-03 1.5383e-02  0.3676    0.7161    
#> d:(Intercept) 3.0520e-01 1.1147e-02 27.3785 < 2.2e-16 ***
#> e:(Intercept) 1.8290e+02 2.6858e+01  6.8098 3.153e-07 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.04562085 (26 degrees of freedom)
#> 
#> Non-normality/heterogeneity adjustment through Box-Cox transformation
#> 
#> Estimated lambda: 0.707 
#> Confidence interval for lambda: [0.439,1.016] 
#> 
```
