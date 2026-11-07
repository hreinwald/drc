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

  a list of non-linear functions to be compared. Models with a log(ED)
  parameterisation (`LL2.*`, `LN.*(loge = TRUE)`) are not allowed (an
  error with instructions is raised); use their dose-scale equivalents
  `LL.*` / `LN.*`.

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
#> LL.4            1.463706 3.057955 6.388640 0.14047306
#> LL.5            1.560325 3.023550 7.729717 0.06816146
#> LN.4            1.489059 3.044614 6.225189 0.12248824
#> W1.4            1.405944 3.088955 5.101078 0.03782468
#> W2.4            1.628230 2.996923 7.806230 0.17886712
#> FPL.4(-1,1)     1.540362 3.038754 7.085943 0.18043369
#> FPL.4(-2,3)     1.507056 3.063616 5.836845 0.08869757
#> FPL.4(-0.5,0.5) 1.531613 3.047967 7.204860 0.18305419
#> 
#>      Estimate
#> e:10 1.530747
#> e:50 3.039441
#> e:90 6.891174

## Model-averaging with Buckland confidence intervals
maED(
  ryegrass.m1,
  list(LL.5(), LN.4(), W1.4(), W2.4()),
  c(10, 50, 90),
  interval = "buckland"
)
#>          ED10     ED50     ED90    Weight
#> LL.4 1.463706 3.057955 6.388640 0.2564245
#> LL.5 1.560325 3.023550 7.729717 0.1244243
#> LN.4 1.489059 3.044614 6.225189 0.2235944
#> W1.4 1.405944 3.088955 5.101078 0.0690465
#> W2.4 1.628230 2.996923 7.806230 0.3265103
#> 
#>      Estimate Std. Error    Lower    Upper
#> e:10 1.531127  0.1960471 1.146882 1.915372
#> e:50 3.032904  0.1943257 2.652033 3.413775
#> e:90 6.892912  1.5548273 3.845507 9.940318

## Model-averaging with Kang confidence intervals
maED(
  ryegrass.m1,
  list(LL.5(), LN.4(), W1.4(), W2.4()),
  c(10, 50, 90),
  interval = "kang"
)
#>          ED10     ED50     ED90    Weight
#> LL.4 1.463706 3.057955 6.388640 0.2564245
#> LL.5 1.560325 3.023550 7.729717 0.1244243
#> LN.4 1.489059 3.044614 6.225189 0.2235944
#> W1.4 1.405944 3.088955 5.101078 0.0690465
#> W2.4 1.628230 2.996923 7.806230 0.3265103
#> 
#>      Estimate Std. Error    Lower    Upper
#> e:10 1.531127  0.1960471 1.155949 1.906305
#> e:50 3.032904  0.1943257 2.631977 3.433831
#> e:90 6.892912  1.5548273 4.241304 9.544520
```
