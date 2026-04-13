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
#> b:(Intercept)   -0.571777    0.187101 -3.0560  0.006001 ** 
#> d:(Intercept)    0.704855    0.025553 27.5845 < 2.2e-16 ***
#> e:(Intercept) 1145.654281  392.624568  2.9179  0.008223 ** 
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
#> e:10  121.8002    59.4612   -1.8561  245.4565
#> e:20  262.9046    72.2694  112.6121  413.1970
#> e:50 1145.6543   392.6246  329.1468 1962.1618

## Model with ED50 as parameter
p.prom.m2<-drm(dryweight~conc, data=P.promelas, fct=LN.3(loge=TRUE))
summary(p.prom.m2)
#> 
#> Model fitted: Log-normal with lower limit at 0 (3 parms)
#> 
#> Parameter estimates:
#> 
#>                Estimate Std. Error t-value   p-value    
#> b:(Intercept) -0.571671   0.187187  -3.054  0.006028 ** 
#> d:(Intercept)  0.704852   0.025555  27.581 < 2.2e-16 ***
#> e:(Intercept)  7.044103   0.343117  20.530 2.221e-15 ***
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
#> e:10  121.79469    0.48838   44.11055  336.29021
#> e:20  262.93031    0.27492  148.43715  465.73480
#> e:50 1146.08015    0.34312  561.46613 2339.41042
```
