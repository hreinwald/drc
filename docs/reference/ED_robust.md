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
  interval = get_ed_interval(mod$fct$name, small_n = TRUE),
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
  [`get_ed_interval()`](https://hreinwald.github.io/drc/reference/get_ed_interval.md).
  Common options include "delta", "tfls", or "buckland".

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
#>     Estimate    stderr     Lower    Upper confint_level confint_method
#>        <num>     <num>     <num>    <num>         <num>         <char>
#> 1:  4.457785  1.674585  1.930237 10.29503          0.95           tfls
#> 2: 35.022556 15.426732 13.125303 93.45151          0.95           tfls
#>           model    EC
#>          <char> <num>
#> 1: BC.4:b-d-e-f    10
#> 2: BC.4:b-d-e-f    50
```
