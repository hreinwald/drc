# Robust Calculation of Model-Averaged Effective Doses

This function serves as a robust wrapper for
[`drc::maED`](https://hreinwald.github.io/drc/reference/maED.md). It
calculates model-averaged effective doses (EDs) for specified response
levels. The key feature is its resilience to errors; it iterates through
each response level individually and handles failures gracefully by
returning `NA` values for that level, rather than terminating the entire
operation.

## Usage

``` r
maED_robust(
  mod,
  fct_ls = NULL,
  respLev = c(10, 20, 50),
  interval = "buckland",
  CI_level = 0.95,
  verbose = FALSE,
  ...
)
```

## Arguments

- mod:

  A model object of class 'drc', which serves as the base model for the
  averaging.

- fct_ls:

  A list of alternative dose-response functions (e.g.,
  [`LL.3()`](https://hreinwald.github.io/drc/reference/LL.3.md),
  [`W1.4()`](https://hreinwald.github.io/drc/reference/W1.4.md)) to be
  used in the model averaging process. The list should be named.

- respLev:

  A numeric vector specifying the response levels (in percentages) for
  which to calculate the EDs (e.g., `c(10, 50)` for EC10 and EC50).

- interval:

  A character string specifying the type of confidence interval to be
  supplied. The default is "buckland". See
  [`drc::maED`](https://hreinwald.github.io/drc/reference/maED.md) for
  other options.

- CI_level:

  A numeric value between 0 and 1 specifying the confidence level for
  the confidence intervals. Default is 0.95.

- verbose:

  A logical value. If `TRUE`, the function will print status messages
  about the calculation progress and any errors encountered for each
  response level. Default is `FALSE`.

- ...:

  Additional arguments to be passed to the underlying
  [`drc::maED`](https://hreinwald.github.io/drc/reference/maED.md)
  function.

## Value

A `data.frame` with one row for each response level specified in
`respLev`. The columns are:

- Estimate:

  The estimated model-averaged effective dose.

- stderr:

  The standard error of the estimate.

- Lower:

  The lower bound of the confidence interval.

- Upper:

  The upper bound of the confidence interval.

- confint_level:

  The confidence level used for the interval.

- confint_method:

  The method used for the confidence interval calculation.

- model:

  A character string listing the models used for averaging.

- EC:

  The response level (as a percentage).

If the calculation for a specific response level fails or results in a
non-positive estimate, the corresponding row will contain `NA` values
for `Estimate`, `stderr`, `Lower`, and `Upper`.

## Details

The function enhances
[`drc::maED`](https://hreinwald.github.io/drc/reference/maED.md) by
introducing a robust calculation loop. It iterates over each element of
`respLev` and calls
[`drc::maED`](https://hreinwald.github.io/drc/reference/maED.md) within
a `tryCatch` block. This approach isolates failures, preventing an error
at one response level (e.g., an EC99 that cannot be estimated) from
halting the calculation of others.

Furthermore, after a successful calculation, the function checks if the
resulting 'Estimate' is positive. If the estimate is `NA`, non-positive,
or if the `tryCatch` block catches an error, the function returns a
structured row of `NA`s for that response level, ensuring a consistent
output format.

## See also

[`maED`](https://hreinwald.github.io/drc/reference/maED.md)

## Author

Hannes Reinwald

## Examples

``` r
data(lettuce)
base_model <- drm(weight ~ conc, data = lettuce, fct = BC.5())
model_list <- list(W2.4 = W2.4())
maED_robust(base_model, fct_ls = model_list, respLev = c(10, 50))
#>     Estimate    stderr      Lower    Upper confint_level confint_method
#>        <num>     <num>      <num>    <num>         <num>         <char>
#> 1:  3.561851  1.610667   0.405001  6.71870          0.95       buckland
#> 2: 11.952400 11.849870 -11.272918 35.17772          0.95       buckland
#>        model    EC
#>       <char> <num>
#> 1: BC.5/W2.4    10
#> 2: BC.5/W2.4    50
```
