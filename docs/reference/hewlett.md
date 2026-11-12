# Hewlett Mixture Model

Provides the Hewlett model for describing the joint action of two
compounds in binary mixture experiments. Used internally by
[`mixture`](https://hreinwald.github.io/drc/reference/mixture.md).

## Usage

``` r
hewlett(
  fixed = c(NA, NA, NA, NA, NA, NA),
  names = c("b", "c", "d", "e", "f", "g"),
  method = c("1", "2", "3", "4"),
  ssfct = NULL,
  eps = 1e-10
)
```

## Arguments

- fixed:

  numeric vector. Specifies which parameters are fixed and at what value
  they are fixed. NAs for parameters that are not fixed.

- names:

  a vector of character strings giving the names of the parameters
  (should not contain ":").

- method:

  character string indicating the self starter function to use.

- ssfct:

  a self starter function to be used (optional).

- eps:

  numeric tolerance for handling zero dose values.

## Value

A list containing the nonlinear model function, the self starter
function, and the parameter names.

## See also

[`mixture`](https://hreinwald.github.io/drc/reference/mixture.md),
[`voelund`](https://hreinwald.github.io/drc/reference/voelund.md)

## Author

Christian Ritz
