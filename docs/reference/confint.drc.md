# Confidence Intervals for Model Parameters

Computes confidence intervals for one or more parameters in a fitted
dose-response model of class `"drc"`. Confidence intervals are
constructed using either a t-distribution (for continuous response
models) or a standard normal distribution (for all other response
types).

## Usage

``` r
# S3 method for class 'drc'
confint(object, parm, level = 0.95, pool = TRUE, ...)
```

## Arguments

- object:

  A fitted model object of class `"drc"`.

- parm:

  A specification of which parameters are to be given confidence
  intervals, either a vector of indices or a vector of parameter name
  strings. If missing, all parameters are considered.

- level:

  The confidence level required. Defaults to `0.95`.

- pool:

  Logical. If `TRUE` (default), curves are pooled. Otherwise they are
  not. This argument only works for models with independently fitted
  curves as specified in
  [`drm()`](https://hreinwald.github.io/drc/reference/drm.md).

- ...:

  Additional arguments for methods. Currently not used.

## Value

A numeric matrix with two columns giving the lower and upper confidence
limits for each parameter. Columns are labelled as \\\frac{(1 -
\text{level})}{2} \times 100\\\\ and \\\left(1 - \frac{(1 -
\text{level})}{2}\right) \times 100\\\\ (by default `2.5 %` and
`97.5 %`).

## See also

- [`drm()`](https://hreinwald.github.io/drc/reference/drm.md) — for
  fitting dose-response models.

- [`confint.basic()`](https://hreinwald.github.io/drc/reference/confint.basic.md)
  — the internal helper used to construct the intervals.

- [`summary.drc()`](https://hreinwald.github.io/drc/reference/summary.drc.md)
  — for a full summary of model coefficients.

## Author

Christian Ritz, Hannes Reinwald

## Examples

``` r
## Fitting a four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

## Confidence intervals for all parameters
confint(ryegrass.m1)
#>                    2.5 %    97.5 %
#> b:(Intercept) 2.01211606 3.9523221
#> c:(Intercept) 0.03878752 0.9240389
#> d:(Intercept) 7.39961398 8.1863026
#> e:(Intercept) 2.67052621 3.4453837

## Confidence interval for a single parameter
confint(ryegrass.m1, "e")
#>                  2.5 %   97.5 %
#> e:(Intercept) 2.670526 3.445384
```
