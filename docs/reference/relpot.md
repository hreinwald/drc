# Relative potency function

Calculates and optionally plots relative potency as a function of the
response level for two curves in a dose-response model, using
[`EDcomp`](https://hreinwald.github.io/drc/reference/EDcomp.md) for the
underlying comparisons.

## Usage

``` r
relpot(
  object,
  plotit = TRUE,
  compMatch = NULL,
  percVec = NULL,
  interval = "none",
  type = c("relative", "absolute"),
  scale = c("original", "percent", "unconstrained"),
  ...
)
```

## Arguments

- object:

  an object of class 'drc'.

- plotit:

  logical. If TRUE (default), a plot of relative potency against
  response level is produced.

- compMatch:

  a numeric vector of length 2 specifying which two curves to compare.

- percVec:

  numeric vector of response levels at which to evaluate relative
  potency. If NULL, a suitable range is determined automatically.

- interval:

  character string specifying confidence interval type. Default is
  "none".

- type:

  character string. Either "relative" (default) or "absolute" response
  levels.

- scale:

  character string. One of "original" (default), "percent", or
  "unconstrained".

- ...:

  additional graphical arguments passed to `plot`.

## Value

An invisible list with components `x`, `y` (relative potency values),
and `percVec`.

## Author

Christian Ritz
