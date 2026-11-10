# Dose-response model selection

Re-fits a dose-response model with alternative dose-response functions
and compares the resulting fits. The original model is always included.

## Usage

``` r
mselect(
  object,
  fctList = NULL,
  nested = FALSE,
  sorted = c("IC", "Res var", "Lack of fit", "no"),
  linreg = FALSE,
  icfct = AIC
)
```

## Arguments

- object:

  A fitted object of class `"drc"`.

- fctList:

  An optional list of dose-response functions to fit to the same data.
  The original model is included even when it is not repeated in this
  list.

- nested:

  Logical; if `TRUE`, report p-values from tests comparing each model
  with the preceding successfully fitted model (starting with `object`).
  These tests are meaningful only when each pair is nested and is
  otherwise comparable.

- sorted:

  Criterion used to order the result. `"IC"` (the default) orders by the
  information criterion in ascending order; `"Res var"` orders by
  residual variance in ascending order; `"Lack of fit"` orders by
  lack-of-fit p-value in descending order; and `"no"` preserves the
  input order. Missing values are placed last.

- linreg:

  Logical; if `TRUE`, also fit linear, quadratic, and cubic polynomial
  regressions. This option is available only for continuous data with a
  single dose-response curve.

- icfct:

  Function used to calculate the information criterion. Defaults to
  [`AIC`](https://rdrr.io/r/stats/AIC.html);
  [`BIC`](https://rdrr.io/r/stats/AIC.html) can also be used. The
  function must return one numeric value for a fitted model.

## Value

A numeric matrix with one row per model (including the original model)
and columns for the log likelihood, information criterion, lack-of-fit
p-value, and, for continuous data, residual variance. If
`nested = TRUE`, a nested-test p-value column is also included unless
`linreg = TRUE`.

## Details

The returned matrix reports the log likelihood, the selected information
criterion, the lack-of-fit p-value, and (for continuous data) the
residual variance. Smaller information-criterion and residual-variance
values indicate a better fit; a larger lack-of-fit p-value indicates
less evidence of lack of fit. A lack-of-fit p-value is `NA` when the fit
type does not provide this test. Log likelihoods should only be compared
for models fit to the same data under comparable likelihoods.

Candidate models are re-fitted to the data stored in `object`. A
candidate that cannot be fitted is retained in the result with `NA`
statistics and a warning. `mselect()` returns a comparison table, not
the fitted candidate objects. With the default `sorted = "IC"`, the
first row is the model with the smallest information criterion.

For nested comparisons, the order of `fctList` determines the sequence
of adjacent comparisons. Failed fits have `NA` statistics; a nested
comparison is omitted if either model in that comparison failed. A test
error or invalid p-value is reported as `NA` with a warning. These tests
should not be used to compare non-nested models.

## Author

Christian Ritz

## Examples

``` r
## Fit a starting model to the ryegrass data
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

## Compare candidate dose-response functions; the smallest AIC is ranked first
candidate_fcts <- list(LL.3 = LL.3(), W1.4 = W1.4(), W2.4 = W2.4())
selection <- mselect(ryegrass.m1, fctList = candidate_fcts)
selection
#>         logLik       IC Lack of fit   Res var
#> W2.4 -15.91352 41.82703   0.9450714 0.2646283
#> LL.4 -16.15514 42.31029   0.8664830 0.2700107
#> W1.4 -17.46720 44.93439   0.4505676 0.3012075
#> LL.3 -18.60413 45.20827   0.3531679 0.3153724

## Identify and refit the AIC-selected model
best_name <- rownames(selection)[1]
best_fit <- if (best_name == ryegrass.m1$fct$name) {
  ryegrass.m1
} else {
  update(ryegrass.m1, fct = candidate_fcts[[best_name]])
}
best_name
#> [1] "W2.4"
best_fit$fct$name
#> [1] "W2.4"

## Add polynomial regressions to the comparison
mselect(ryegrass.m1, fctList = candidate_fcts, linreg = TRUE)
#>          logLik        IC Lack of fit   Res var
#> W2.4  -15.91352  41.82703   0.9450714 0.2646283
#> LL.4  -16.15514  42.31029   0.8664830 0.2700107
#> W1.4  -17.46720  44.93439   0.4505676 0.3012075
#> LL.3  -18.60413  45.20827   0.3531679 0.3153724
#> Cubic -25.53428  61.06856          NA 0.5899609
#> Quad  -35.11558  78.23116          NA 1.2485122
#> Lin   -50.47554 106.95109          NA 4.2863247

## Compare adjacent models in the returned order
mselect(ryegrass.m1, fctList = list(LL.5()), nested = TRUE)
#>         logLik       IC Lack of fit   Res var Nested F test
#> LL.4 -16.15514 42.31029   0.8664830 0.2700107            NA
#> LL.5 -15.87828 43.75656   0.8538476 0.2777393     0.5134602
```
