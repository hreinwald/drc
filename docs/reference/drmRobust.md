# Robust estimation functions for drm

Robust estimation functions for drm

## Usage

``` r
drmRobust(robust, fctCall, lenData, lenPar, callEnv = parent.frame())
```

## Arguments

- robust:

  character string, the rho function (see
  [`drm`](https://hreinwald.github.io/drc/reference/drm.md)).

- fctCall:

  the (matched) call of
  [`drm()`](https://hreinwald.github.io/drc/reference/drm.md).

- lenData:

  number of observations.

- lenPar:

  number of estimated parameters.

- callEnv:

  environment in which `fctCall` is evaluated for the auxiliary LTS fit
  (scale estimate of `"trimmed"`, `"tukey"` and `"winsor"`);
  [`drm()`](https://hreinwald.github.io/drc/reference/drm.md) passes the
  environment it was called from.
