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

The parameters are: `b`, the slope (positive for a decreasing curve);
`c`, the lower asymptote; `d`, the upper asymptote; and `e`, the ED50
(or log(ED50) if `loge = TRUE`).

If `loge = TRUE`,
[`ED`](https://hreinwald.github.io/drc/reference/ED.md) returns
**log(ED)**; use `interval = "fls"` to obtain EDs and confidence limits
on the dose scale (`"fls"` is only valid for such models). Model
averaging with
[`maED`](https://hreinwald.github.io/drc/reference/maED.md) requires
`loge = FALSE`.

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

## Examples

``` r
## Fit with drm() and plot
m1 <- drm(rootl ~ conc, data = ryegrass, fct = lnormal())
plot(m1, type = "all", main = "Four-parameter log-normal model")


## Effect of the lower asymptote c: c free, fixed at 0 or fixed at 3
m.c0 <- drm(rootl ~ conc, data = ryegrass, fct = LN.4(fixed = c(NA, 0, NA, NA)))
m.c3 <- drm(rootl ~ conc, data = ryegrass, fct = LN.4(fixed = c(NA, 3, NA, NA)))
plot(m1, type = "all", main = "Lower asymptote c")
plot(m.c0, type = "none", add = TRUE, col = "darkgreen")
plot(m.c3, type = "none", add = TRUE, col = "red")
legend("topright", c("c free", "c = 0", "c = 3"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")


## Effect of the upper asymptote d: d free, fixed at 6 or fixed at 9
m.d6 <- drm(rootl ~ conc, data = ryegrass, fct = LN.4(fixed = c(NA, NA, 6, NA)))
m.d9 <- drm(rootl ~ conc, data = ryegrass, fct = LN.4(fixed = c(NA, NA, 9, NA)))
plot(m1, type = "all", main = "Upper asymptote d")
plot(m.d6, type = "none", add = TRUE, col = "darkgreen")
plot(m.d9, type = "none", add = TRUE, col = "red")
legend("topright", c("d free", "d = 6", "d = 9"), col = c("black", "darkgreen", "red"),
  lty = 1, bty = "n")


## log(ED50) as parameter: ED() returns log(ED), 'fls' back-transforms
m2 <- drm(rootl ~ conc, data = ryegrass, fct = LN.4(loge = TRUE))
ED(m2, c(10, 50), interval = "fls")
#> 
#> Estimated effective doses
#> 
#>      Estimate Std. Error    Lower    Upper
#> e:10 1.489060   0.120623 1.157810 1.915081
#> e:50 3.044615   0.060995 2.680870 3.457713
```
