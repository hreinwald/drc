# Fractional polynomial-logistic dose-response model

Model function for specifying dose-response models that are a
combination of a logistic model and an appropriate class of fractional
polynomials.

## Usage

``` r
fplogistic(
  p1,
  p2,
  fixed = c(NA, NA, NA, NA),
  names = c("b", "c", "d", "e"),
  method = c("1", "2", "3", "4"),
  ssfct = NULL,
  fctName,
  fctText
)
```

## Arguments

- p1:

  numeric denoting the negative power of log(dose+1) in the fractional
  polynomial.

- p2:

  numeric denoting the positive power of log(dose+1) in the fractional
  polynomial.

- fixed:

  numeric vector. Specifies which parameters are fixed and at what value
  they are fixed. NAs for parameters that are not fixed.

- names:

  a vector of character strings giving the names of the parameters
  (should not contain ":"). The order of the parameters is: b, c, d, e.

- method:

  character string indicating the self starter function to use.

- ssfct:

  a self starter function to be used.

- fctName:

  optional character string used internally by convenience functions.

- fctText:

  optional character string used internally by convenience functions.

## Value

A list containing the nonlinear function, the self starter function and
the parameter names.

## Details

The fractional polynomial dose-response models introduced by Namata et
al. (2008) are implemented using the logistic model as base.

## References

Namata, Harriet and Aerts, Marc and Faes, Christel and Teunis, Peter
(2008) Model Averaging in Microbial Risk Assessment Using Fractional
Polynomials, *Risk Analysis* **28**, 891–905.

## See also

[`FPL.4`](https://hreinwald.github.io/drc/reference/FPL.4.md),
[`maED`](https://hreinwald.github.io/drc/reference/maED.md),
[`drm`](https://hreinwald.github.io/drc/reference/drm.md)

## Author

Christian Ritz
