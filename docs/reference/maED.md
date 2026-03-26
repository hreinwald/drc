# Estimation of ED values using model-averaging

Estimates and confidence intervals for ED values are estimated using
model-averaging.

## Usage

``` r
maED(
  object,
  fctList = NULL,
  respLev = c(10, 20, 50),
  interval = c("none", "buckland", "kang"),
  linreg = FALSE,
  clevel = NULL,
  level = 0.95,
  type = c("relative", "absolute"),
  display = TRUE,
  na.rm = FALSE,
  extended = FALSE
)
```

## Arguments

- object:

  an object of class `drc`.

- fctList:

  a list of non-linear functions to be compared.

- respLev:

  a numeric vector containing the response levels.

- interval:

  character string specifying the type of confidence intervals to be
  supplied. The default is `"none"`. The choices `"buckland"` and
  `"kang"` are explained in the Details section.

- linreg:

  logical indicating whether or not additionally a simple linear
  regression model should be fitted.

- clevel:

  character string specifying the curve id in case estimates for a
  specific curve or compound are requested. By default estimates are
  shown for all curves.

- level:

  numeric. The confidence level. Must be a single value strictly between
  0 and 1. The default is `0.95`.

- type:

  character string. Whether the specified response levels are absolute
  or relative (default).

- display:

  logical. If `TRUE` results are displayed. Otherwise they are not
  (useful in simulations).

- na.rm:

  logical indicating whether or not `NA` values occurring during model
  fitting should be excluded from subsequent calculations.

- extended:

  logical specifying whether or not an extended output (including fit
  summaries) should be returned.

## Value

If `extended = FALSE`, a matrix with two or more columns containing the
model-averaged estimates and the corresponding estimated standard errors
and, optionally, lower and upper confidence limits. If
`extended = TRUE`, a list with components:

- estimates:

  Matrix of model-averaged ED estimates and intervals.

- fits:

  Matrix of per-model ED estimates and AIC-based weights.

## Details

Model-averaging of individual estimates is carried out as described by
Buckland *et al.* (1997) and Kang *et al.* (2000) using AIC-based
weights. The two approaches differ w.r.t. the calculation of confidence
intervals: Buckland *et al.* (1997) provide an approximate variance
formula under the assumption of perfectly correlated estimates (so,
confidence intervals will tend to be too wide). Kang *et al.* (2000) use
the model weights to calculate confidence limits as weighted means of
the confidence limits for the individual fits.

## References

Buckland, S. T. and Burnham, K. P. and Augustin, N. H. (1997) Model
Selection: An Integral Part of Inference, *Biometrics* **53**, 603–618.

Kang, Seung-Ho and Kodell, Ralph L. and Chen, James J. (2000)
Incorporating Model Uncertainties along with Data Uncertainties in
Microbial Risk Assessment, *Regulatory Toxicology and Pharmacology*
**32**, 68–72.

## See also

The function
[`mselect`](https://hreinwald.github.io/drc/reference/mselect.md)
provides a summary of fit statistics for several models fitted to the
same data.

## Author

Christian Ritz, Hannes Reinwald

## Examples

``` r
## Fitting an example dose-response model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

## Model-averaging with default settings (no confidence intervals)
maED(
  ryegrass.m1,
  list(LL.5(), LN.4(), W1.4(), W2.4(), FPL.4(-1, 1), FPL.4(-2, 3), FPL.4(-0.5, 0.5)),
  c(10, 50, 90)
)
#>                     ED10     ED50     ED90     Weight
#> LL.4            1.463706 3.057955 6.388640 0.14047308
#> LL.5            1.560325 3.023549 7.729713 0.06816147
#> LN.4            1.489188 3.044673 6.224889 0.12248817
#> W1.4            1.405979 3.088964 5.101022 0.03782468
#> W2.4            1.628278 2.996913 7.805803 0.17886712
#> FPL.4(-1,1)     1.540346 3.038790 7.086271 0.18043370
#> FPL.4(-2,3)     1.507055 3.063612 5.836831 0.08869758
#> FPL.4(-0.5,0.5) 1.531613 3.047967 7.204860 0.18305421
#> 
#>      Estimate
#> e:10 1.530770
#> e:50 3.039453
#> e:90 6.891117

## Model-averaging with Buckland confidence intervals
maED(
  ryegrass.m1,
  list(LL.5(), LN.4(), W1.4(), W2.4()),
  c(10, 50, 90),
  interval = "buckland"
)
#>          ED10     ED50     ED90     Weight
#> LL.4 1.463706 3.057955 6.388640 0.25642453
#> LL.5 1.560325 3.023549 7.729713 0.12442435
#> LN.4 1.489188 3.044673 6.224889 0.22359424
#> W1.4 1.405979 3.088964 5.101022 0.06904651
#> W2.4 1.628278 2.996913 7.805803 0.32651037
#> 
#>      Estimate Std. Error    Lower    Upper
#> e:10 1.531174  0.1977800 1.143532 1.918816
#> e:50 3.032914  0.1945023 2.651697 3.414132
#> e:90 6.892701  1.5391238 3.876074 9.909329

## Model-averaging with Kang confidence intervals
maED(
  ryegrass.m1,
  list(LL.5(), LN.4(), W1.4(), W2.4()),
  c(10, 50, 90),
  interval = "kang"
)
#>          ED10     ED50     ED90     Weight
#> LL.4 1.463706 3.057955 6.388640 0.25642453
#> LL.5 1.560325 3.023549 7.729713 0.12442435
#> LN.4 1.489188 3.044673 6.224889 0.22359424
#> W1.4 1.405979 3.088964 5.101022 0.06904651
#> W2.4 1.628278 2.996913 7.805803 0.32651037
#> 
#>      Estimate    Lower    Upper
#> e:10 1.531174 1.155988 1.906360
#> e:50 3.032914 2.631988 3.433841
#> e:90 6.892701 4.241165 9.544238
```
