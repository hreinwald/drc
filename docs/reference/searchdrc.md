# Search through a range of initial parameter values to obtain convergence

`searchdrc` provides a facility for searching through a range of initial
values for a single parameter in order to obtain convergence of the
non-linear estimation procedure used in dose-response curve fitting.

## Usage

``` r
searchdrc(object, which, range, len = 50, verbose = FALSE)
```

## Arguments

- object:

  an object of class `'drc'`, which must have valid `$start` and
  `$parNames` fields populated. This is typically an object from a model
  that failed to converge but was still constructed with initial
  parameter values.

- which:

  a character string containing the parameter name **without** the curve
  suffix (e.g., `"b"` not `"b:1"`). Must exactly match one of the
  parameter names in the model object.

- range:

  a numeric vector of exactly length 2 specifying the interval endpoints
  `c(lower, upper)` for the search range. The two endpoints must be
  different.

- len:

  a positive integer (minimum 2). The maximum number of evenly spaced
  starting values to try within `range`. The search stops early as soon
  as convergence is achieved, so the actual number of attempts may be
  less than `len`. Defaults to `50`.

- verbose:

  logical. If `TRUE`, prints progress messages indicating which starting
  value is currently being tried. Defaults to `FALSE`.

## Value

If convergence is achieved, returns the fitted model object of class
`'drc'`, corresponding to the **first** starting value in the search
grid that led to a successful fit. If no starting value leads to
convergence, the function throws an error.

## Details

The function iterates through at most `len` evenly spaced values within
the specified `range`, using each as a starting value for the chosen
parameter. The search stops as soon as the first successful model fit is
found. You would need to identify the parameter which is most likely to
cause problems for the estimation procedure.

Parameter names should be provided **without** the curve suffix. For
example, use `"b"` rather than `"b:1"`. The function internally matches
the parameter using the pattern `"^<which>:"` against the full parameter
names stored in the model object.

## See also

[`drm`](https://hreinwald.github.io/drc/reference/drm.md) for the main
model fitting function,
[`drmc`](https://hreinwald.github.io/drc/reference/drmc.md) for control
arguments, [`update`](https://rdrr.io/r/stats/update.html) for the
update method used internally.

## Author

Christian Ritz, Hannes Reinwald.

## Examples

``` r
if (FALSE) { # \dontrun{
library(drc)

# Fit an initial model (which may fail to converge)
myModel <- drm(response ~ dose, data = myData, fct = LL.4())

# Search over a range of starting values for the slope parameter "b"
myModelFixed <- searchdrc(myModel, which = "b", range = c(-5, 5), len = 100)

# With progress messages enabled
myModelFixed <- searchdrc(myModel, which = "b", range = c(-5, 5),
                          len = 100, verbose = TRUE)
} # }
```
