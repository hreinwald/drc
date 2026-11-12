# Fitting binary mixture models

`mixture` fits a concentration addition, Hewlett or Voelund model to
data from binary mixture toxicity experiments.

## Usage

``` r
mixture(
  object,
  model = c("CA", "Hewlett", "Voelund"),
  start,
  startm,
  control = drmc()
)
```

## Arguments

- object:

  object of class 'drc' corresponding to the model with freely varying
  EC50 values.

- model:

  character string. It can be "CA", "Hewlett" or "Voelund".

- start:

  optional numeric vector supplying starting values for all parameters
  in the mixture model.

- startm:

  optional numeric vector supplying the lambda parameter in the Hewlett
  model or the eta parameters (two parameters) in the Voelund model.

- control:

  list of arguments controlling constrained optimisation (zero as
  boundary), maximum number of iteration in the optimisation, relative
  tolerance in the optimisation, warnings issued during the
  optimisation.

## Value

An object of class 'drc' with a few additional components.

## Details

The function is a wrapper to
[`drm`](https://hreinwald.github.io/drc/reference/drm.md), implementing
the models described in Soerensen et al. (2007). See the paper for a
discussion of the merits of the different models.

Currently only the log-logistic models are available. Application of
Box-Cox transformation is not yet available.

## References

Ritz, C. and Streibig, J. C. (2014) From additivity to synergism - A
modelling perspective *Synergy*, **1**, 22–29.

## Author

Christian Ritz
