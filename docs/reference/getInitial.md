# Showing starting values used

Returns the starting values of the model parameters used when fitting a
dose-response model.

## Usage

``` r
getInitial(object)
```

## Arguments

- object:

  object of class 'drc'.

## Value

A vector of starting values for the model parameters used to initialize
the estimation procedure.

## Note

This function is masking the standard function in the stats package.

## Author

Christian Ritz
