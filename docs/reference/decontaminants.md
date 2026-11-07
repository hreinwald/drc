# Performance of decontaminants used in the culturing of a micro-organism

The two decontaminants 1-hexadecylpyridium chloride and oxalic acid were
used. Additionally there was a control group (coded as concentration 0
and only included under oxalic acid).

## Usage

``` r
data("decontaminants")
```

## Format

A data frame with 128 observations on the following 3 variables.

- `conc`:

  a numeric vector of percentage weight per volume

- `count`:

  a numeric vector of numbers of M. bovis colonies at stationarity

- `group`:

  a factor with levels `hpc` and `oxalic` of the decontaminants used

## Details

These data examplify Wadley's problem: counts where the maximum number
is not known. The data were analyzed by Trajstman (1989) using a
three-parameter logistic model and then re-analyzed by Morgan and Smith
(1992) using a three-parameter Weibull type II model. In both cases the
authors adjusted for overdispersion (in different ways).

It seems that Morgan and Smith (1992) fitted separate models for the two
decontaminants and using the control group for both model fits. In the
example below a joint model is fitted where the control group is used
once to determine a shared upper limit at concentration 0.

## Source

Trajstman, A. C. (1989) Indices for Comparing Decontaminants when Data
Come from Dose-Response Survival and Contamination Experiments, *Applied
Statistics*, **38**, 481–494.

## References

Morgan, B. J. T. and Smith, D. M. (1992) A Note on Wadley's Problem with
Overdispersion, *Applied Statistics*, **41**, 349–354.

## Examples

``` r
library(drc)

## Wadley's problem using a three-parameter log-logistic model
decon.LL.3.1 <- drm(count~conc, group, data = decontaminants, fct = LL.3(), 
type = "Poisson", pmodels = list(~group, ~1, ~group))

summary(decon.LL.3.1)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) with lower limit at 0 (3 parms)
#> 
#> Parameter estimates:
#> 
#>                Estimate Std. Error t-value   p-value    
#> b:(Intercept)  1.271768   0.085526 14.8700 < 2.2e-16 ***
#> b:groupoxalic -0.750855   0.096728 -7.7626 8.269e-15 ***
#> d:(Intercept) 48.956059   1.162693 42.1058 < 2.2e-16 ***
#> e:(Intercept)  0.161702   0.011731 13.7847 < 2.2e-16 ***
#> e:groupoxalic  0.346366   0.109175  3.1726  0.001511 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

plot(decon.LL.3.1)



## Same model fit in another parameterization (no intercepts)
decon.LL.3.2 <- drm(count~conc, group, data = decontaminants, fct=LL.3(), 
type = "Poisson", pmodels = list(~group-1, ~1, ~group-1))

summary(decon.LL.3.2)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) with lower limit at 0 (3 parms)
#> 
#> Parameter estimates:
#> 
#>                Estimate Std. Error t-value   p-value    
#> b:grouphpc     1.271673   0.085517 14.8703 < 2.2e-16 ***
#> b:groupoxalic  0.521034   0.059805  8.7123 < 2.2e-16 ***
#> d:(Intercept) 48.955627   1.162671 42.1062 < 2.2e-16 ***
#> e:grouphpc     0.161729   0.011734 13.7831 < 2.2e-16 ***
#> e:groupoxalic  0.507902   0.113186  4.4873 7.213e-06 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```
