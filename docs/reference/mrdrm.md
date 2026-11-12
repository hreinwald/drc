# Model-robust dose-response modelling

Combines a parametric dose-response curve fitted by
[`drm()`](https://hreinwald.github.io/drc/reference/drm.md) with a
non-parametric local-linear regression (`loess`) fit. The fitted
response is a weighted average of the two curves. When `lambda` contains
more than one value, the function selects the weight that minimizes a
prediction-error criterion; this reduces reliance on a single parametric
curve shape when that shape may be misspecified. The result is intended
as a complementary, model-robust fit, not as a replacement for the
interpretable parameters of the parametric
[`drm()`](https://hreinwald.github.io/drc/reference/drm.md) fit.

## Usage

``` r
mrdrm(
  object1,
  object2,
  lambda = (0:10)/10,
  criterion = c("gcv", "lcv"),
  critFct = c("ls", "ll"),
  ls.weights = c("nonpar", "ad hoc", "none", "par", "response"),
  fixedEnd = FALSE,
  unitScale = FALSE
)
```

## Arguments

- object1:

  A single-curve model fitted by
  [`drm()`](https://hreinwald.github.io/drc/reference/drm.md), with
  `type = "continuous"` or `type = "binomial"`.

- object2:

  An optional degree-1 `loess` fit to the same response and dose data.
  If omitted, one is fitted using the formula and data stored in
  `object1`.

- lambda:

  Numeric candidate mixing weights between 0 and 1, inclusive. The
  combined fit is `(1 - lambda) * object1 + lambda * object2`, so 0 is
  purely parametric and 1 is purely local-linear. A single value fixes
  the weight and skips selection; multiple values are compared using
  `criterion`.

- criterion:

  Prediction-error criterion used to select `lambda`: `"gcv"` for
  generalized cross-validation or `"lcv"` for leave-one-dose- level-out
  cross-validation. With a single `lambda`, no criterion is evaluated.

- critFct:

  Criterion function: `"ls"` for least squares or `"ll"` for binomial
  log-likelihood. Continuous fits always use least squares.

- ls.weights:

  Weighting scheme for the least-squares cross-validation criterion on
  binomial data: `"nonpar"` (weights from local-linear predictions),
  `"par"` (weights from parametric predictions), `"response"` (weights
  from observed proportions), `"ad hoc"` (a boundary-adjusted scheme),
  or `"none"`. This argument does not affect continuous fits.

- fixedEnd:

  If `TRUE`, use observed mean responses at the smallest and largest
  dose levels for the local-linear leave-one-dose-out predictions. This
  tends to favor the non-parametric component at the endpoints.

- unitScale:

  If `TRUE`, refit the local-linear regression using dose values
  rescaled to the unit interval. Any supplied `object2` is replaced by
  this internally fitted regression. This option cannot be combined with
  leave-one-dose-out selection.

## Value

An object of class `"mrdrc"`, containing the selected `lambda`, the
combined fitted values, goodness-of-fit summaries, and both component
models. Use [`predict()`](https://rdrr.io/r/stats/predict.html) for
predictions and
[`ED()`](https://hreinwald.github.io/drc/reference/ED.md) for
effective-dose estimates.

## Details

Use this function when the parametric model is scientifically reasonable
but its exact shape is uncertain, or as a sensitivity analysis alongside
the parametric fit. Effective doses for the resulting fit are obtained
by inverse regression with
[`ED()`](https://hreinwald.github.io/drc/reference/ED.md).

Only single-curve continuous and binomial
[`drm()`](https://hreinwald.github.io/drc/reference/drm.md) fits are
supported. Continuous fits always use least squares. For binomial fits,
least squares or a binomial log-likelihood can be used. The
cross-validation criterion leaves out one dose level at a time (all
replicate observations at that dose), rather than one individual
observation.

## References

Nottingham, S. F. and Birch, C. P. D. (2000). A semi-parametric approach
to analysing dose-response data. *Environmental and Ecological
Statistics*, 7, 85-100.

## See also

[`drm()`](https://hreinwald.github.io/drc/reference/drm.md),
[`ED()`](https://hreinwald.github.io/drc/reference/ED.md),
[`predict()`](https://rdrr.io/r/stats/predict.html)

## Examples

``` r
# \donttest{
data(ryegrass)
fit <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
robustFit <- mrdrm(fit)
robustFit$lambda
#> [1] 0
ED(robustFit, c(10, 50, 90), display = FALSE)
# }
```
