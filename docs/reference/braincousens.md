# The Brain-Cousens hormesis models

`braincousens` provides a very general way of specifying Brain-Cousens'
modified log-logistic model for describing hormesis, under various
constraints on the parameters.

## Usage

``` r
braincousens(
  fixed = c(NA, NA, NA, NA, NA),
  names = c("b", "c", "d", "e", "f"),
  method = c("1", "2", "3", "4"),
  ssfct = NULL,
  fctName,
  fctText
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

## Value

A list containing the non-linear function, the self starter function,
the parameter names and additional model specific objects.

## Details

The Brain-Cousens model is given by the expression \$\$f(x) = c +
\frac{d-c+fx}{1+\exp(b(\log(x)-\log(e)))}\$\$ which is a five-parameter
model.

## References

Brain, P. and Cousens, R. (1989) An equation to describe dose responses
where there is stimulation of growth at low doses, *Weed Research*,
**29**, 93–96.

## See also

[`BC.4`](https://hreinwald.github.io/drc/reference/BC.4.md),
[`BC.5`](https://hreinwald.github.io/drc/reference/BC.5.md),
[`drm`](https://hreinwald.github.io/drc/reference/drm.md)

## Author

Christian Ritz
