# Searching through a range of initial parameter values to obtain convergence

`searchdrc` provides a facility for searching through a range of
parameter values (one-dimensional) in order to obtain convergence of the
estimation procedure.

## Usage

``` r
searchdrc(object, which, range, len = 50)
```

## Arguments

- object:

  an object of class 'drc'. The object can be from a model that could
  not be fitted.

- which:

  a character string containing the parameter name.

- range:

  a numeric vector of length 2 specifying the interval endpoints for the
  range.

- len:

  numeric. The number of points in the interval.

## Value

An object of class 'drc'.

## Details

The function goes through the range with increments such that in total
at most `len` sets of parameter values are used as initial values for
the estimation procedure. You would need to identify the parameter which
is most likely to cause problems for the estimation procedure.

## Author

Christian Ritz
