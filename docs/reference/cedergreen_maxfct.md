# Find the Dose and Response at Maximum Hormesis

This function finds the dose that elicits the maximum hormetic
(stimulatory) response for the Cedergreen-Ritz model and the response
value at that dose.

## Usage

``` r
cedergreen_maxfct(
  all_params,
  alpha,
  lower = 1e-06,
  upper = 1000,
  .optimize_fn = stats::optimize
)
```

## Arguments

- all_params:

  A named list of all model parameters (b, c, d, e, f).

- alpha:

  The hormesis alpha shape parameter.

- lower:

  The lower bound of the dose interval to search for the maximum.

- upper:

  The upper bound of the dose interval to search for the maximum.

## Value

A numeric vector containing two values: the dose at the maximum
response, and the maximum response value itself. Returns `c(NA, NA)` on
failure.

## Author

Hannes Reinwald
