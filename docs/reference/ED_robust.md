# Robust Calculation of Effective Doses (ED)

This function serves as a robust wrapper for
[`drc::ED`](https://hreinwald.github.io/drc/reference/ED.md). It
calculates effective doses (EDs) for multiple specified response levels.
Its primary feature is the ability to gracefully handle cases where an
ED value is not mathematically estimable from the model (e.g., the
requested response is outside the model's asymptotes). Instead of
throwing an error, it returns a row of `NA` values for that specific
response level, ensuring the overall analysis can proceed.

## Usage

``` r
ED_robust(
  mod,
  respLev = c(10, 20, 50),
  interval = get_ed_interval(mod$fct$name, small_n = TRUE, log_ed =
    .isLogEDModel(mod$fct)),
  CI_level = 0.95,
  verbose = FALSE,
  ...
)
```

## Arguments

- mod:

  An object of class 'drc', representing the fitted dose-response model.

- respLev:

  A numeric vector specifying the response levels for which to calculate
  ED values (e.g., `c(10, 50)` for ED10 and ED50).

- interval:

  A character string specifying the method for calculating confidence
  intervals. Defaults to the output of
  [`get_ed_interval()`](https://hreinwald.github.io/drc/reference/get_ed_interval.md):
  "tfls" for the log-logistic, log-normal, Brain-Cousens and Cedergreen
  families, "delta" for Weibull models, and "fls" for models returning
  log(ED) (`LL2.*`, `LN.*(loge = TRUE)`), for which estimates and limits
  are back-transformed and the standard error is the delta-method SE of
  the dose. "fls" is only accepted for such models (an error is raised
  otherwise).

- CI_level:

  A numeric value between 0 and 1 indicating the confidence level for
  the intervals (e.g., 0.95 for a 95% CI).

- verbose:

  A logical value. If `TRUE`, the function will print status messages
  about the calculation progress and any errors encountered for each
  response level. Default is `FALSE`.

- ...:

  Additional arguments to be passed directly to
  [`drc::ED`](https://hreinwald.github.io/drc/reference/ED.md).

## Value

A `data.table` where each row corresponds to a requested response level.
The table includes the ED estimate, standard error, confidence interval
(Lower, Upper), and metadata about the calculation (confidence level,
method, model name, and EC level). Rows for non-estimable EDs are
populated with `NA`.

## Author

Hannes Reinwald

## Examples

``` r
data(lettuce)
m <- drm(weight ~ conc, data = lettuce, fct = BC.4())
ED_robust(m, respLev = c(10, 50), CI_level = 0.95)
#>     Estimate   stderr     Lower    Upper confint_level confint_method
#>        <num>    <num>     <num>    <num>         <num>         <char>
#> 1:  4.458442  1.67505  1.930311 10.29767          0.95           tfls
#> 2: 35.045444 15.44581 13.126369 93.56610          0.95           tfls
#>           model    EC
#>          <char> <num>
#> 1: BC.4:b-d-e-f    10
#> 2: BC.4:b-d-e-f    50
```
