# Creating isobolograms

`isobole` displays isobole based on EC/ED50 estimates from a
log-logistic model. Additionally isoboles determined by the
concentration addition model, Hewlett's model and Voelund's model can be
added to the plot.

## Usage

``` r
isobole(
  object1,
  object2,
  exchange = 1,
  cifactor = 2,
  ename = "e",
  xaxis = "100",
  xlab,
  ylab,
  xlim,
  ylim,
  ...
)
```

## Arguments

- object1:

  object of class 'drc' where EC/ED50 parameters vary freely.

- object2:

  object of class 'drc' where EC/ED50 parameters vary according to
  Hewlett's model.

- exchange:

  numeric. The exchange rate between the two substances.

- cifactor:

  numeric. The factor to be used in the confidence intervals. Default is
  2, but 1 has been used in publications.

- ename:

  character string. The name of the EC/ED50 variable.

- xaxis:

  character string. Is the mixture "0:100" or "100:0" on the x axis?

- xlab:

  an optional label for the x axis.

- ylab:

  an optional label for the y axis.

- xlim:

  a numeric vector of length two, containing the lower and upper limit
  for the x axis.

- ylim:

  a numeric vector of length two, containing the lower and upper limit
  for the y axis.

- ...:

  Additional graphical parameters.

## Value

No value is returned. Only used for the side effect: the isobologram
shown.

## Details

The model fits to be supplied as first and optionally second argument
are obtained using
[`mixture`](https://hreinwald.github.io/drc/reference/mixture.md) and
[`drm`](https://hreinwald.github.io/drc/reference/drm.md).

## References

Ritz, C. and Streibig, J. C. (2014) From additivity to synergism - A
modelling perspective *Synergy*, **1**, 22–29.

## Author

Christian Ritz
