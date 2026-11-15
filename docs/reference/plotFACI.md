# Plot combination index as a function of fraction affected

Visualizes the combination index from
[`CIcompX`](https://hreinwald.github.io/drc/reference/CIcompX.md) as a
function of the fraction affected.

## Usage

``` r
plotFACI(
  effList,
  indAxis = c("ED", "EF"),
  caRef = TRUE,
  showPoints = FALSE,
  add = FALSE,
  ylim,
  ...
)
```

## Arguments

- effList:

  a list as returned by
  [`CIcompX`](https://hreinwald.github.io/drc/reference/CIcompX.md).

- indAxis:

  character string. Either "ED" for effective doses or "EF" for effects.

- caRef:

  logical. If TRUE (default), a reference line for concentration
  addition is drawn.

- showPoints:

  logical. If TRUE, estimated combination indices are plotted as points.

- add:

  logical. If TRUE, the plot is added to an existing plot.

- ylim:

  numeric vector of length 2 giving the range for the y axis.

- ...:

  additional graphical arguments.

## Value

Invisibly returns the plot matrix of combination index values.

## See also

[`CIcompX`](https://hreinwald.github.io/drc/reference/CIcompX.md),
[`CIcomp`](https://hreinwald.github.io/drc/reference/CIcomp.md)

## Author

Christian Ritz and Ismael Rodea-Palomares
