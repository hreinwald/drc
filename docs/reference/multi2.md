# Multistage Dose-Response Model with Quadratic Terms

A five-parameter multistage dose-response model useful for describing
more complex dose-response patterns.

## Usage

``` r
multi2(
  fixed = c(NA, NA, NA, NA, NA),
  names = c("b1", "b2", "b3", "c", "d"),
  ssfct = NULL,
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

- ssfct:

  a self starter function to be used.

- fctName:

  optional character string used internally by convenience functions.

- fctText:

  optional character string used internally by convenience functions.

## Value

A list containing the nonlinear function, the self starter function, and
the parameter names.

## Details

The multistage model function with quadratic terms is:

\$\$f(x) = c + (d-c)\exp(-b1 - b2 x - b3 x^2)\$\$

where x denotes the dose or the logarithm-transformed dose.

The parameters are: `b1`, `b2`, `b3`, the polynomial coefficients (`b1`
sets the response at dose 0 relative to the range, `b2` and `b3` the
linear and quadratic decline); `c`, the lower asymptote; and `d`, the
upper bound of the scale. No effective doses
([`ED`](https://hreinwald.github.io/drc/reference/ED.md)) are available
for this model.

## References

Wheeler, M. W., Bailer, A. J. (2009) Comparing model averaging with
other model selection strategies for benchmark dose estimation,
*Environmental and Ecological Statistics*, **16**, 37–51.

## See also

[`gammadr`](https://hreinwald.github.io/drc/reference/gammadr.md),
[`drm`](https://hreinwald.github.io/drc/reference/drm.md)

## Author

Christian Ritz

## Examples

``` r
## Fit with drm() and plot
m1 <- drm(rootl ~ conc, data = ryegrass, fct = multi2())
plot(m1, type = "all", main = "Multistage model with quadratic terms")

```
