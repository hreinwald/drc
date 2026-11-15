# Select Appropriate Confidence Interval Method for a drc Model

This function determines the recommended confidence interval calculation
method ('type' argument in drc::ED) based on the model family of a 'drc'
object.

## Usage

``` r
get_ed_interval(
  model,
  small_n = TRUE,
  fls_pattern = "^LL|^LN|^BC|^CRS",
  verbose = FALSE,
  log_ed = FALSE
)
```

## Arguments

- model:

  A drc model object or a character string specifying the model name
  (e.g., "LL.4").

- small_n:

  A logical value. If TRUE, "tfls" (confidence limits computed on the
  log scale of the ED and back-transformed) is used for models matched
  by `fls_pattern`; if FALSE the delta method ("delta") is used. Models
  with log(ED) parameterisation (`LL2.*`) always use "fls", the only
  valid choice for them. Defaults to TRUE.

- fls_pattern:

  A regular expression character string. This pattern is used to
  identify model families for which the "fls" or "tfls" method is
  appropriate. The default covers standard log-logistic, log-normal,
  Brain-Cousens, and Cedergreen-Ritz-Streibig models.

- verbose:

  A logical value. If TRUE, a message is printed when the function
  resorts to its default choice because the model type was not
  explicitly matched. Defaults to TRUE.

- log_ed:

  Logical; `TRUE` if the model returns ED on the log scale (`LL2.*`,
  `LN.*(loge = TRUE)`). Detected automatically for `drc` objects and for
  `LL2.*` names.

## Value

A character string: "tfls", "fls", or "delta", representing the
recommended interval type for use in
[`drc::ED()`](https://hreinwald.github.io/drc/reference/ED.md).

## Author

Hannes Reinwald

## Examples

``` r
ryegrass_model <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
drc:::get_ed_interval(ryegrass_model)
#> [1] "tfls"
drc:::get_ed_interval("LL.4")
#> [1] "tfls"
drc:::get_ed_interval("W1.4")
#> [1] "delta"
```
