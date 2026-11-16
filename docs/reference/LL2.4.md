# Four-Parameter Log-Logistic Model with log(ED50) as Parameter

A four-parameter log-logistic model where the ED50 is parameterised on
the log scale. The asymmetry parameter `f` is fixed at 1. The estimated
parameters are the slope `b`, the lower limit `c`, the upper limit `d`,
and the log(ED50) `e`.

## Usage

``` r
LL2.4(fixed = c(NA, NA, NA, NA), names = c("b", "c", "d", "e"), ...)
```

## Arguments

- fixed:

  numeric vector of length 4. Specifies which parameters are fixed and
  at what value. Use `NA` for parameters that should be estimated.

- names:

  character vector of length 4 giving the names of the parameters `b`,
  `c`, `d`, and `e`.

- ...:

  additional arguments passed to
  [`llogistic2`](https://hreinwald.github.io/drc/reference/llogistic2.md).

## Value

A list of class `"llogistic"` with the nonlinear function, self-starter,
and related components.

## See also

[`llogistic2`](https://hreinwald.github.io/drc/reference/llogistic2.md),
[`LL2.3`](https://hreinwald.github.io/drc/reference/LL2.3.md),
[`LL2.5`](https://hreinwald.github.io/drc/reference/LL2.5.md)

## Examples

``` r
## Several curves fitted with drm() and plotted; ED on the dose scale
spinach.m1 <- drm(SLOPE ~ DOSE, CURVE, data = spinach, fct = LL2.4())
plot(spinach.m1, main = "LL2.4, one curve per CURVE level")

ED(spinach.m1, 50, interval = "fls")
#> 
#> Estimated effective doses
#> 
#>        Estimate Std. Error    Lower    Upper
#> e:1:50 1.794886   0.266433 1.056753 3.048600
#> e:2:50 0.945678   0.263895 0.559592 1.598141
#> e:3:50 1.372715   0.329661 0.712721 2.643875
#> e:4:50 0.197327   0.051642 0.178072 0.218665
#> e:5:50 0.210173   0.065310 0.184579 0.239315
```
