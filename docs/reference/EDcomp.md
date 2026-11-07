# Comparison of relative potencies between dose-response curves

Relative potencies (also called selectivity indices) for arbitrary doses
are compared between fitted dose-response curves.

## Usage

``` r
EDcomp(
  object,
  percVec,
  percMat = NULL,
  compMatch = NULL,
  od = FALSE,
  vcov. = vcov,
  reverse = FALSE,
  interval = c("none", "delta", "fieller", "fls"),
  level = ifelse(!(interval == "none"), 0.95, NULL),
  reference = c("control", "upper"),
  type = c("relative", "absolute"),
  display = TRUE,
  pool = TRUE,
  logBase = NULL,
  multcomp = FALSE,
  ...
)
```

## Arguments

- object:

  an object of class 'drc'.

- percVec:

  a numeric vector of dosage values.

- percMat:

  a matrix with 2 columns providing the pairs of indices of `percVec` to
  be compared. By default all pairs are compared.

- compMatch:

  an optional character vector of names of assays to be compared. If not
  specified all comparisons are supplied.

- od:

  logical. If TRUE adjustment for over-dispersion is used. This argument
  only makes a difference for binomial data.

- vcov.:

  function providing the variance-covariance matrix.
  [`vcov`](https://rdrr.io/r/stats/vcov.html) is the default, but
  `sandwich` is also an option (for obtaining robust standard errors).

- reverse:

  logical. If TRUE the order of comparison of two curves is reversed.

- interval:

  character string specifying the type of confidence intervals to be
  supplied. The default is `"none"`. Use `"delta"` for asymptotics-based
  confidence intervals, `"fieller"` for confidence intervals based on
  Fieller's theorem, or `"fls"` for confidence intervals
  back-transformed from logarithm scale.

- level:

  numeric. The level for the confidence intervals. Default is 0.95.

- reference:

  character string. Is the upper limit or the control level the
  reference?

- type:

  character string specifying whether absolute or relative response
  levels are supplied.

- display:

  logical. If TRUE results are displayed. Otherwise they are not (useful
  in simulations).

- pool:

  logical. If TRUE curves are pooled. Otherwise they are not. This
  argument only works for models with independently fitted curves as
  specified in
  [`drm`](https://hreinwald.github.io/drc/reference/drm.md).

- logBase:

  numeric. The base of the logarithm in case logarithm transformed dose
  values are used.

- multcomp:

  logical to switch on output for use with the package multcomp. Default
  is FALSE.

- ...:

  additional arguments passed to the function doing the calculations.

## Value

An invisible matrix containing the estimates and the corresponding
estimated standard errors and possibly lower and upper confidence
limits. Or, alternatively, a list with elements that may be plugged
directly into `parm` in the package multcomp (when `multcomp` is TRUE).

## Details

Fieller's theorem is incorporated using the formulas provided by Kotz
and Johnson (1983) and Finney (1978).

For objects of class 'braincousens' or 'mlogistic' the additional
argument may be the 'upper' argument or the 'interval' argument
specifying limits for the bisection method.

## See also

[`ED.drc`](https://hreinwald.github.io/drc/reference/ED.drc.md) for
calculating effective doses.

## Author

Christian Ritz

## Examples

``` r
spinach.LL.4 <- drm(SLOPE~DOSE, CURVE, data = spinach, fct = LL.4())

EDcomp(spinach.LL.4, c(50, 50))
#> 
#> Estimated ratios of effect doses
#> 
#>             Estimate Std. Error    t-value    p-value
#> 1/2:50/50  1.8983289  0.7118572  1.2619511  0.2104182
#> 1/3:50/50  1.3072639  0.5541673  0.5544606  0.5807197
#> 1/4:50/50  9.0964443  2.4687287  3.2796006  0.0015079
#> 1/5:50/50  8.5152036  2.3364916  3.2164479  0.0018366
#> 2/3:50/50  0.6886393  0.2908179 -1.0706378  0.2873630
#> 2/4:50/50  4.7918167  1.2884228  2.9429911  0.0041895
#> 2/5:50/50  4.4856313  1.2196574  2.8578774  0.0053619
#> 3/4:50/50  6.9583840  2.3222736  2.5657545  0.0120489
#> 3/5:50/50  6.5137602  2.1897834  2.5179478  0.0136790
#> 4/5:50/50  0.9361024  0.0781404 -0.8177281  0.4157988
EDcomp(spinach.LL.4, c(10, 50))
#> 
#> Estimated ratios of effect doses
#> 
#>              Estimate  Std. Error     t-value     p-value
#> 1/2:10/50  2.7644e-02  1.5249e-02 -6.3765e+01  7.3857e-74
#> 1/3:10/50  1.9037e-02  1.1155e-02 -8.7938e+01  1.5164e-85
#> 1/4:10/50  1.3247e-01  6.4532e-02 -1.3444e+01  4.7646e-23
#> 1/5:10/50  1.2400e-01  6.0616e-02 -1.4452e+01  6.4884e-25
#> 2/3:10/50  4.4296e-02  3.7848e-02 -2.5251e+01  1.4383e-41
#> 2/4:10/50  3.0823e-01  2.4348e-01 -2.8412e+00  5.6247e-03
#> 2/5:10/50  2.8853e-01  2.2822e-01 -3.1174e+00  2.4893e-03
#> 3/4:10/50  2.7741e-01  1.5001e-01 -4.8170e+00  6.2896e-06
#> 3/5:10/50  2.5969e-01  1.4081e-01 -5.2574e+00  1.0720e-06
#> 4/5:10/50  2.8449e-01  3.7366e-02 -1.9149e+01  7.0655e-33
EDcomp(spinach.LL.4, c(10, 50), reverse = TRUE)
#> 
#> Estimated ratios of effect doses
#> 
#>             Estimate Std. Error    t-value    p-value
#> 2/1:50/10 3.6174e+01 1.9954e+01 1.7627e+00 8.1543e-02
#> 3/1:50/10 5.2530e+01 3.0781e+01 1.6741e+00 9.7796e-02
#> 4/1:50/10 7.5491e+00 3.6776e+00 1.7808e+00 7.8518e-02
#> 5/1:50/10 8.0644e+00 3.9422e+00 1.7920e+00 7.6689e-02
#> 3/2:50/10 2.2576e+01 1.9289e+01 1.1185e+00 2.6650e-01
#> 4/2:50/10 3.2444e+00 2.5629e+00 8.7573e-01 3.8365e-01
#> 5/2:50/10 3.4658e+00 2.7414e+00 8.9948e-01 3.7094e-01
#> 4/3:50/10 3.6047e+00 1.9492e+00 1.3363e+00 1.8502e-01
#> 5/3:50/10 3.8508e+00 2.0881e+00 1.3653e+00 1.7577e-01
#> 5/4:50/10 3.5151e+00 4.6169e-01 5.4476e+00 4.8854e-07
```
