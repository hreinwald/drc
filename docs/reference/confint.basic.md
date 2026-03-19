# Basic Confidence Interval Calculation

An internal helper function that constructs a confidence interval matrix
from a matrix of parameter estimates and their standard errors. A
t-distribution quantile is used for continuous response models; a
standard normal quantile is used for all other response types (binomial,
event, Poisson, negbin1, negbin2).

## Usage

``` r
# S3 method for class 'basic'
confint(estMat, level, intType, dfres, formatting = TRUE)
```

## Arguments

- estMat:

  A numeric matrix with two columns: the first column contains parameter
  estimates and the second column contains their standard errors.

- level:

  The confidence level required (e.g., `0.95` for 95% intervals).

- intType:

  A character string specifying the response type of the model. One of
  `"binomial"`, `"continuous"`, `"event"`, `"Poisson"`, `"negbin1"`, or
  `"negbin2"`. Determines whether a normal or t-distribution quantile is
  used. For `"continuous"` models a t-distribution with `dfres` degrees
  of freedom is used; all other types use the standard normal.

- dfres:

  The residual degrees of freedom. Only used when
  `intType = "continuous"`.

- formatting:

  Logical. If `TRUE` (default), row and column names are added to the
  returned matrix.

## Value

A numeric matrix with two columns giving the lower and upper confidence
limits for each parameter.

## See also

[`confint.drc()`](https://hreinwald.github.io/drc/reference/confint.drc.md)
— the user-facing function that calls this helper.
