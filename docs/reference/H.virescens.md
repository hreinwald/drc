# Mortality of tobacco budworms

For three days, moths of the tobacco budworm (*Heliothis virescens*)
were exposed to doses of the pyrethroid trans-cypermethrin.

## Usage

``` r
data(H.virescens)
```

## Format

A data frame with 12 observations on the following 4 variables.

- `dose`:

  a numeric vector of dose values (\\\mu g\\)

- `numdead`:

  a numeric vector of dead or knocked-down moths

- `total`:

  a numeric vector of total number of moths

- `sex`:

  a factor with levels `F` `M` denoting a grouping according to sex

## Details

In Venables and Ripley (2002), these data are analysed using a logistic
regression with base-2 logarithm of dose as explanatory variable.

## Source

Venables, W. N. and Ripley, B. D (2002) *Modern Applied Statistics with
S*, New York: Springer (fourth edition).

## Examples

``` r
library(drc)

## Fitting dose-response model (log-logistic with common slope)
Hv.m1 <- drm(numdead/total~dose, sex, weights = total, data = H.virescens, fct = LL.2(), 
pmodels = list(~ 1, ~ sex - 1), type = "binomial")
summary(Hv.m1)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) with lower limit at 0 and upper limit at 1 (2 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) -1.53537    0.18911 -8.1189 4.573e-16 ***
#> e:sexF         9.60556    1.52990  6.2786 3.417e-10 ***
#> e:sexM         4.69001    0.73465  6.3840 1.725e-10 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

## Fitting the same model as in Venables and Riply (2002)
Hv.m2 <- glm(cbind(numdead, total-numdead) ~ sex + I(log2(dose)) - 1, data = H.virescens, 
family = binomial)

## Comparing the fits
logLik(Hv.m1)
#> 'log Lik.' -18.43373 (df=3)
logLik(Hv.m2)
#> 'log Lik.' -18.43373 (df=3)

## Estimated ED values (matching those given in MASS)
ED(Hv.m1, c(25, 50, 75))
#> 
#> Estimated effective doses
#> 
#>        Estimate Std. Error
#> e:F:25  4.69645    0.81353
#> e:F:50  9.60556    1.52990
#> e:F:75 19.64607    3.74120
#> e:M:25  2.29309    0.41882
#> e:M:50  4.69001    0.73465
#> e:M:75  9.59239    1.69565
```
