# Self-starter for the Cedergreen-Ritz-Streibig Dose-Response Model

A self-starting function for the Cedergreen-Ritz-Streibig model, used to
find initial parameter estimates for non-linear regression (e.g., with
`nls` or `drc`).

## Usage

``` r
cedergreen.ssf(
  method = c("loglinear", "anke", "method3", "normolle"),
  fixed,
  alpha,
  useFixed = FALSE
)
```

## Arguments

- method:

  A character string specifying the method for estimating initial 'b'
  and 'e' parameters. Using descriptive names is preferred.

- fixed:

  A numeric vector of fixed parameter values, with `NA` for parameters
  that need to be estimated. The required order is `c(b, c, d, e, f)`.

- alpha:

  A numeric value for the alpha parameter, which is treated as a known
  constant during the estimation of the other initial parameters.

- useFixed:

  A logical value. If `TRUE`, the function will use the non-NA values
  provided in the `fixed` argument as fixed parameters and only estimate
  the others.

## Value

A numeric vector of initial parameter estimates for the model parameters
that were not specified as `fixed`.

## Details

This function is a closure that returns another function. The returned
function takes a data frame and calculates initial values for the model
parameters (b, c, d, e, f). This self-starter relies on several helper
functions (e.g., `findcd`, `findbe1`, `findbe2`, `findbe3`) which must
be available in the calling environment.
