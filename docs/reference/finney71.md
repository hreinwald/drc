# Example from Finney (1971)

For each of six concentrations of an insecticide the number of insects
affected (out of the total number of insects) was recorded.

## Usage

``` r
data(finney71)
```

## Format

A data frame with 6 observations on the following 3 variables.

- `dose`:

  a numeric vector

- `total`:

  a numeric vector

- `affected`:

  a numeric vector

## Source

Finney, D. J. (1971) *Probit Analysis*, Cambridge: Cambridge University
Press.

## Examples

``` r
library(drc)

## Model with ED50 as a parameter
finney71.m1 <- drm(affected/total ~ dose, weights = total,
data = finney71, fct = LL.2(), type = "binomial")

summary(finney71.m1)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) with lower limit at 0 and upper limit at 1 (2 parms)
#> 
#> Parameter estimates:
#> 
#>               Estimate Std. Error t-value   p-value    
#> b:(Intercept) -3.10363    0.38773 -8.0047 1.154e-15 ***
#> e:(Intercept)  4.82890    0.24958 19.3485 < 2.2e-16 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
plot(finney71.m1, broken = TRUE, bp = 0.1, lwd = 2)


ED(finney71.m1, c(10, 20, 50), interval = "delta", reference = "control")
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error   Lower   Upper
#> e:10  2.37896    0.25164 1.88576 2.87217
#> e:20  3.08932    0.24372 2.61163 3.56700
#> e:50  4.82890    0.24958 4.33974 5.31806

## Model fitted with 'glm'
#fitl.glm <- glm(cbind(affected, total-affected) ~ log(dose),
#family=binomial(link = logit), data=finney71[finney71$dose != 0, ])
#summary(fitl.glm)  # p-value almost agree for the b parameter
#
#xp <- dose.p(fitl.glm, p=c(0.50, 0.90, 0.95))  # from MASS
#xp.ci <- xp + attr(xp, "SE") %*% matrix(qnorm(1 - 0.05/2)*c(-1,1), nrow=1)
#zp.est <- exp(cbind(xp.ci[,1],xp,xp.ci[,2]))
#dimnames(zp.est)[[2]] <- c("zp.lcl","zp","zp.ucl")
#zp.est  # not far from above results with 'ED'

## Model with log(ED50) as a parameter
finney71.m2 <- drm(affected/total ~ dose, weights = total,
data = finney71, fct = LL2.2(), type = "binomial")

## Confidence intervals based on back-transformation
##  complete agreement with results based on 'glm'
ED(finney71.m2, c(10, 20, 50), interval = "fls", reference = "control")
#> 
#> Estimated effective doses
#> 
#>      Estimate  Lower  Upper
#> e:10   2.3789 1.9335 2.9270
#> e:20   3.0893 2.6467 3.6059
#> e:50   4.8289 4.3637 5.3437
```
