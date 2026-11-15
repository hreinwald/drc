# Interactive dose-response modelling

Internal helper used by
[`drm`](https://hreinwald.github.io/drc/reference/drm.md) when
`separate = TRUE`: the model is fitted to each curve separately and the
fits are combined into one `drc` object.

## Usage

``` r
idrm(
  x,
  y,
  curveid,
  weights,
  fct,
  type,
  control,
  robust = "mean",
  bcVal = NULL,
  bcAdd = 0,
  lowerl = NULL,
  upperl = NULL,
  logDose = NULL,
  start = NULL
)
```

## Arguments

- x:

  dose values.

- y:

  response values.

- curveid:

  curve identifiers.

- weights:

  weights.

- fct:

  model function (list).

- type:

  response type.

- control:

  list from [`drmc`](https://hreinwald.github.io/drc/reference/drmc.md).

- robust, bcVal, bcAdd, lowerl, upperl, logDose, start:

  arguments of [`drm`](https://hreinwald.github.io/drc/reference/drm.md)
  that are forwarded to the per-curve fits.
