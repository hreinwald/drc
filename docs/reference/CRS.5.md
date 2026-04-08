# Wrapper for 5-parameter Cedergreen-Ritz-Streibig Model

A convenience wrapper for the
[`drc::cedergreen`](https://hreinwald.github.io/drc/reference/cedergreen.md)
function, preset for a 5-parameter model. It provides flexible handling
for the alpha parameter.

## Usage

``` r
CRS.5(
  names = c("b", "c", "d", "e", "f"),
  fixed = c(NA, NA, NA, NA, NA),
  alpha_type = "a",
  fctName = NULL,
  fctText = NULL,
  ...
)
```

## Arguments

- names:

  A character vector of length 5 specifying the names of the model
  parameters. Default is `c("b", "c", "d", "e", "f")`.

- fixed:

  A numeric vector of length 5. Use `NA` for parameters to be estimated
  and a numeric value for parameters to be fixed. Default is all `NA`.

- alpha_type:

  A character or a numeric value. Can be one of 'a' (alpha=1), 'b'
  (alpha=0.5), 'c' (alpha=0.25), or a specific numeric value for alpha.

- fctName:

  An optional character string to name the model function. If `NULL`
  (the default), a name is generated automatically.

- fctText:

  An optional character string describing the model. If `NULL` (the
  default), a description is generated automatically.

- ...:

  Additional arguments to be passed to
  [`drc::cedergreen`](https://hreinwald.github.io/drc/reference/cedergreen.md),
  such as `data`.

## Value

A `drc` model object of class `cedergreen`. If the underlying
[`drc::cedergreen`](https://hreinwald.github.io/drc/reference/cedergreen.md)
call fails, it issues a warning and returns `NULL`.

## Details

This function simplifies the creation of a 5-parameter
Cedergreen-Ritz-Streibig model by setting sensible defaults for the
parameter names. It allows the alpha parameter to be specified either by
a predefined character shortcut ('a', 'b', 'c') or by a direct numeric
value.

By default the function runs with `alpha=1`, which corresponds to the
`CRS.4a` model. Setting `alpha=0.5` corresponds to the `CRS.4b` model,
and `alpha=0.25` corresponds to the `CRS.4c` model.

By default, all parameters are set to be estimated (i.e., `fixed` is all
`NA`), but users can specify any parameters to be held constant during
estimation. The self-starter function is automatically generated based
on the specified method and fixed parameters, ensuring that initial
values are appropriately calculated for the model fitting process.

The function automatically generates a model name (`fctName`) and
description (`fctText`) unless they are explicitly provided by the user.

## Author

Hannes Reinwald

## Examples

``` r
# Create a CRS.5 model specification
crs_model_a <- CRS.5()

# Fix the lower limit to 0 and use a custom numeric alpha
crs_model_custom <- CRS.5(
  fixed = c(NA, 0, NA, NA, NA), alpha_type = 0.75
)
```
