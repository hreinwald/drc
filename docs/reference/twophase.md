# Two-Phase Dose-Response Model

A seven-parameter dose-response model combining two log-logistic
components, useful for describing more complex dose-response patterns.

## Usage

``` r
twophase(
  fixed = c(NA, NA, NA, NA, NA, NA, NA),
  names = c("b1", "c1", "d1", "e1", "b2", "d2", "e2"),
  fctName,
  fctText
)
```

## Arguments

- fixed:

  numeric vector specifying which parameters are fixed and at what value
  they are fixed. NAs are used for parameters that are not fixed.

- names:

  a vector of character strings giving the names of the parameters
  (should not contain ":"). The default is reasonable.

- fctName:

  optional character string used internally by convenience functions.

- fctText:

  optional character string used internally by convenience functions.

## Value

A list containing the nonlinear function, the self starter function, and
the parameter names.

## Details

Following Groot *et al* (1996) the two-phase model function is:

\$\$f(x) = c + \frac{d1-c}{1+\exp(b1(\log(x)-\log(e1)))} +
\frac{d2}{1+\exp(b2(\log(x)-\log(e2)))}\$\$

For each of the two phases, the parameters have the same interpretation
as in the ordinary log-logistic model.

## References

Groot, J. C. J., Cone, J. W., Williams, B. A., Debersaques, F. M. A.,
Lantinga, E. A. (1996) Multiphasic analysis of gas production kinetics
for in vitro fermentation of ruminant feeds, *Animal Feed Science
Technology*, **64**, 77–89.

## See also

The basic component in the two-phase model is the log-logistic model
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md).

## Author

Christian Ritz
