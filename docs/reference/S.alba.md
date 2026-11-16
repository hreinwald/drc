# Potency of two herbicides

Data are from an experiment, comparing the potency of the two herbicides
glyphosate and bentazone in white mustard *Sinapis alba*.

## Usage

``` r
data(S.alba)
```

## Format

A data frame with 68 observations on the following 3 variables.

- `Dose`:

  a numeric vector containing the dose in g/ha.

- `Herbicide`:

  a factor with levels `Bentazone` `Glyphosate` (the two herbicides
  applied).

- `DryMatter`:

  a numeric vector containing the response (dry matter in g/pot).

## Details

The lower and upper limits for the two herbicides can be assumed
identical, whereas slopes and ED50 values are different (in the
log-logistic model).

## Source

Christensen, M. G. and Teicher, H. B., and Streibig, J. C. (2003)
Linking fluorescence induction curve and biomass in herbicide screening,
*Pest Management Science*, **59**, 1303–1310.

## Examples

``` r
library(drc)

## Fitting a log-logistic model with
##  common lower and upper limits
S.alba.LL.4.1 <- drm(DryMatter~Dose, Herbicide, data=S.alba, fct = LL.4(),
pmodels=data.frame(Herbicide,1,1,Herbicide)) 
summary(S.alba.LL.4.1)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) (4 parms)
#> 
#> Parameter estimates:
#> 
#>                Estimate Std. Error t-value   p-value    
#> b:Bentazone    5.058358   1.045148  4.8398 8.988e-06 ***
#> b:Glyphosate   2.385448   0.493324  4.8355 9.133e-06 ***
#> c:(Intercept)  0.716452   0.089247  8.0278 3.542e-11 ***
#> d:(Intercept)  3.854375   0.076249 50.5500 < 2.2e-16 ***
#> e:Bentazone   28.638392   2.043019 14.0177 < 2.2e-16 ***
#> e:Glyphosate  66.933054   5.979355 11.1940 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.370514 (62 degrees of freedom)

## Applying the optimal transform-both-sides Box-Cox transformation
## (using the initial model fit)  
S.alba.LL.4.2 <- boxcox(S.alba.LL.4.1, method = "anova") 

summary(S.alba.LL.4.2)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) (4 parms)
#> 
#> Parameter estimates:
#> 
#>                Estimate Std. Error t-value   p-value    
#> b:Bentazone    4.838699   0.927280  5.2182 2.217e-06 ***
#> b:Glyphosate   1.944310   0.236471  8.2222 1.630e-11 ***
#> c:(Intercept)  0.682591   0.028768 23.7271 < 2.2e-16 ***
#> d:(Intercept)  3.862606   0.106186 36.3760 < 2.2e-16 ***
#> e:Bentazone   28.396254   1.874642 15.1476 < 2.2e-16 ***
#> e:Glyphosate  65.573487   5.618957 11.6700 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1558947 (62 degrees of freedom)
#> 
#> Non-normality/heterogeneity adjustment through Box-Cox transformation
#> 
#> Estimated lambda: 0.101 
#> Confidence interval for lambda: [-0.126, 0.331] 
#> 

## Plotting fitted regression curves together with the data
plot(S.alba.LL.4.2)
```
