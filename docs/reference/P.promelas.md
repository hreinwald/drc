# Effect of sodium pentachlorophenate on growth of fathead minnow

Fathead minnows (*Pimephales promelas*) were exposed to sodium
pentachlorophenate concentrations ranging from 32 to 512 micro g/L in a
7-day larval survival and growth test. The average dry weight was
measured.

## Usage

``` r
data(P.promelas)
```

## Format

A data frame with 24 observations on the following 2 variables.

- `conc`:

  a numeric vector of sodium pentachlorophenate concentrations (micro
  g/L).

- `dryweight`:

  a numeric vector dry weights (mg)

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

## Model with ED50 on log scale as parameter
p.prom.m1<-drm(dryweight~conc, data=P.promelas, fct=LN.3())

plot(fitted(p.prom.m1), residuals(p.prom.m1))


plot(p.prom.m1, type="all", broken=TRUE, xlim=c(0,1000))

summary(p.prom.m1)
#> 
#> Model fitted: Log-normal with lower limit at 0 (3 parms)
#> 
#> Parameter estimates:
#> 
#>                  Estimate  Std. Error t-value   p-value    
#> b:(Intercept)   -0.571521    0.187068 -3.0551  0.006012 ** 
#> d:(Intercept)    0.704867    0.025554 27.5831 < 2.2e-16 ***
#> e:(Intercept) 1146.311329  393.304953  2.9146  0.008287 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.05467719 (21 degrees of freedom)
ED(p.prom.m1, c(10,20,50), interval="delta")
#> 
#> Estimated effective doses
#> 
#>       Estimate Std. Error     Lower     Upper
#> e:10  121.7478    59.4563   -1.8983  245.3939
#> e:20  262.8820    72.2762  112.5754  413.1886
#> e:50 1146.3113   393.3050  328.3889 1964.2338

## Model with ED50 as parameter
p.prom.m2<-drm(dryweight~conc, data=P.promelas, fct=LN.3(loge=TRUE))
summary(p.prom.m2)
#> 
#> Model fitted: Log-normal with lower limit at 0 (3 parms)
#> 
#> Parameter estimates:
#> 
#>                Estimate Std. Error t-value   p-value    
#> b:(Intercept) -0.571411   0.187041  -3.055  0.006014 ** 
#> d:(Intercept)  0.704874   0.025555  27.583 < 2.2e-16 ***
#> e:(Intercept)  7.044535   0.343240  20.524 2.228e-15 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.05467719 (21 degrees of freedom)
ED(p.prom.m2, c(10,20,50), interval="fls")
#> 
#> Estimated effective doses
#> 
#>        Estimate Std. Error      Lower      Upper
#> e:10  121.72354    0.48842   44.08134  336.12004
#> e:20  262.86841    0.27496  148.38988  465.66383
#> e:50 1146.57534    0.34324  561.56575 2341.01708
```
