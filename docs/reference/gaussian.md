# Normal (Gaussian) biphasic dose-response model

Model function for fitting symmetric or skewed bell-shaped/biphasic
dose-response patterns using the Gaussian (normal distribution) model.

## Usage

``` r
gaussian(
  fixed = c(NA, NA, NA, NA, NA),
  names = c("b", "c", "d", "e", "f"),
  method = c("1", "2", "3", "4"),
  ssfct = NULL,
  fctName,
  fctText,
  loge = FALSE
)
```

## Arguments

- fixed:

  numeric vector. Specifies which parameters are fixed and at what value
  they are fixed. NAs for parameters that are not fixed.

- names:

  a vector of character strings giving the names of the parameters
  (should not contain ":"). The order of the parameters is: b, c, d, e,
  f.

- method:

  character string indicating the self starter function to use.

- ssfct:

  a self starter function to be used.

- fctName:

  optional character string used internally by convenience functions.

- fctText:

  optional character string used internally by convenience functions.

- loge:

  logical indicating whether or not e or log(e) should be a parameter in
  the model. By default e is a model parameter.

## Value

The value returned is a list containing the nonlinear function, the self
starter function and the parameter names.

## See also

[`lgaussian`](https://hreinwald.github.io/drc/reference/lgaussian.md),
[`drm`](https://hreinwald.github.io/drc/reference/drm.md)

## Author

Christian Ritz
