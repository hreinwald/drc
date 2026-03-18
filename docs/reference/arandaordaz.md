# Asymptotic Regression Model

The base function for the asymptotic regression model, providing the
mean function and self starter for a three-parameter model.

## Usage

``` r
arandaordaz(fixed = c(NA, NA, NA), names = c("a", "b", "c"), fctName, fctText)
```

## Arguments

- fixed:

  numeric vector. Specifies which parameters are fixed and at what value
  they are fixed. Use `NA` for parameters that are not fixed. Must be of
  length 3.

- names:

  character vector of length 3 giving the names of the parameters
  (should not contain ":").

- fctName:

  optional character string used internally by convenience functions.
  Defaults to `"arandaordaz"` if not provided.

- fctText:

  optional character string used internally by convenience functions.
  Defaults to `"Asymptotic regression"` if not provided.

## Value

A list of class `drcMean` with the following components:

- fct:

  The mean function taking arguments `dose` and `parm`.

- ssfct:

  Self-starter function for generating initial parameter estimates from
  data.

- names:

  Character vector of non-fixed parameter names.

- deriv1:

  Reserved first derivative slot (currently `NULL`).

- deriv2:

  Reserved second derivative slot (currently `NULL`).

- derivx:

  Reserved derivative-with-respect-to-x slot (currently `NULL`).

- edfct:

  Function for calculating effective dose (ED) values and their
  derivatives.

- inversion:

  Inverse mean function for back-calculating dose from response.

- name:

  Character string identifying the model function name.

- text:

  Character string with a human-readable model description.

- noParm:

  Integer giving the number of non-fixed parameters.

## Details

The asymptotic regression model is a three-parameter model with mean
function:

\$\$f(x) = c + (d-c)(1-\exp(-x/e))\$\$

The parameter \\c\\ is the lower limit (at \\x=0\\), \\d\\ is the upper
limit, and \\e\>0\\ determines the steepness of the increase.

## See also

[`AR.2`](https://hreinwald.github.io/drc/reference/AR.2.md),
[`AR.3`](https://hreinwald.github.io/drc/reference/AR.3.md),
[`EXD.2`](https://hreinwald.github.io/drc/reference/EXD.2.md),
[`EXD.3`](https://hreinwald.github.io/drc/reference/EXD.3.md)

## Author

Christian Ritz, Hannes Reinwald
