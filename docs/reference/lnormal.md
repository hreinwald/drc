# Log-normal dose-response model

`lnormal` provides a general framework for specifying the mean function
of the decreasing or increasing log-normal dose-response model.

## Usage

``` r
lnormal(
  fixed = c(NA, NA, NA, NA),
  names = c("b", "c", "d", "e"),
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

  vector of character strings giving the names of the parameters (should
  not contain ":"). The order of the parameters is: b, c, d, e.

- method:

  character string indicating the self starter function to use.

- ssfct:

  a self starter function to be used.

- fctName:

  optional character string used internally by convenience functions.

- fctText:

  optional character string used internally by convenience functions.

- loge:

  logical indicating whether or not ED50 or log(ED50) should be a
  parameter in the model. By default ED50 is a model parameter.

## Value

A list containing the non-linear function, the self starter function and
the parameter names.

## Details

For the case where log(ED50) is a parameter in the model, the mean
function is: \$\$f(x) = c + (d-c)(\Phi(b(\log(x)-e)))\$\$ and in case
ED50 is a parameter: \$\$f(x) = c + (d-c)(\Phi(b(\log(x)-\log(e))))\$\$

For \\c=0\\ and \\d=1\\, the model reduces to the classic probit model.

## References

Finney, D. J. (1971) *Probit analysis*, London: Cambridge University
Press.

Bruce, R. D. and Versteeg, D. J. (1992) A statistical procedure for
modeling continuous toxicity data, *Environ. Toxicol. Chem.*, **11**,
1485–1494.

## See also

The log-logistic model
[`llogistic`](https://hreinwald.github.io/drc/reference/llogistic.md) is
very similar to the log-normal model.

## Author

Christian Ritz
