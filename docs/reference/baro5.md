# The Baroreflex Five-Parameter Dose-Response Model

`baro5` provides the five-parameter baroreflex model function, allowing
specification under various parameter constraints. The model
accommodates asymmetric dose-response curves.

## Usage

``` r
baro5(
  fixed = c(NA, NA, NA, NA, NA),
  names = c("b1", "b2", "c", "d", "e"),
  method = c("1", "2", "3", "4"),
  ssfct = NULL
)
```

## Arguments

- fixed:

  numeric vector. Specifies which parameters are fixed and at what value
  they are fixed. NAs for parameters that are not fixed.

- names:

  a vector of character strings giving the names of the parameters
  (should not contain ":"). The order is: b1, b2, c, d, e.

- method:

  character string indicating the self starter function to use.

- ssfct:

  a self starter function to be used.

## Value

A list containing the nonlinear model function, the self starter
function, and the parameter names.

## Details

The five-parameter function is given by:

\$\$y = c + \frac{d-c}{1+f\exp(b1(\log(x)-\log(e))) +
(1-f)\exp(b2(\log(x)-\log(e)))}\$\$

\$\$f = 1/(1 + \exp((2b1 b2/\|b1+b2\|)(\log(x)-\log(e))))\$\$

If the difference between b1 and b2 is nonzero, the function is
asymmetric.

## References

Ricketts, J. H. and Head, G. A. (1999) A five-parameter logistic
equation for investigating asymmetry of curvature in baroreflex studies.
*Am. J. Physiol. (Regulatory Integrative Comp. Physiol. 46)*, **277**,
441–454.

## Author

Christian Ritz
