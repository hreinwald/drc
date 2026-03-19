# Universal Response Surface Approach (URSA) for Drug Interaction

URSA provides a parametric approach for modelling the joint action of
several agents. The model allows quantification of synergistic effects
through a single parameter. The model function is defined implicitly
through an appropriate equation.

## Usage

``` r
ursa(
  fixed = rep(NA, 7),
  names = c("b1", "b2", "c", "d", "e1", "e2", "f"),
  ssfct = NULL
)
```

## Arguments

- fixed:

  numeric vector. Specifies which parameters are fixed and at what value
  they are fixed. NAs for parameters that are not fixed.

- names:

  a vector of character strings giving the names of the parameters. The
  default is reasonable.

- ssfct:

  a self starter function to be used (optional).

## Value

A list containing the nonlinear function, the self starter function, and
the parameter names.

## References

Greco, W. R. and Park H. S. and Rustum, Y. M. (1990) Application of a
New Approach for the Quantitation of Drug Synergism to the Combination
of cis-Diamminedichloroplatinum and 1-beta-D-Arabinofuranosylcytosine,
*Cancer Research*, **50**, 5318–5327.

Greco, W. R. Bravo, G. and Parsons, J. C. (1995) The Search for Synergy:
A Critical Review from a Response Surface Perspective, *Pharmacological
Reviews*, **47**, Issue 2, 331–385.

## See also

Other models for fitting mixture data:
[`mixture`](https://hreinwald.github.io/drc/reference/mixture.md).

## Author

Christian Ritz after an idea by Hugo Ceulemans.

## Examples

``` r
d1 <- c(0, 0, 0, 0, 0, 0, 0, 0, 2, 5, 10, 20, 50, 2, 2, 2,
  2, 2, 5, 5, 5, 5, 5, 10, 10, 10, 10, 10, 20, 20, 20, 20,
  20, 50, 50, 50, 50, 50)
d2 <- c(0, 0, 0, 0.2, 0.5, 1, 2, 5, 0, 0, 0, 0, 0, 0.2,
  0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5, 0.2,
  0.5, 1, 2, 5, 0.2, 0.5, 1, 2, 5)
effect <- c(106, 99.2, 115, 79.2, 70.1, 49, 21, 3.83, 74.2,
  71.5, 48.1, 30.9, 16.3, 76.3, 48.8, 44.5, 15.5, 3.21,
  56.7, 47.5, 26.8, 16.9, 3.25, 46.7, 35.6, 21.5, 11.1,
  2.94, 24.8, 21.6, 17.3, 7.78, 1.84, 13.6, 11.1, 6.43,
  3.34, 0.89)
greco <- data.frame(d1, d2, effect)
greco.m1 <- drm(effect ~ d1 + d2, data = greco,
  fct = ursa(fixed = c(NA, NA, 0, NA, NA, NA, NA)))
summary(greco.m1)
#> 
#> Model fitted: URSA (6 parms)
#> 
#> Parameter estimates:
#> 
#>                  Estimate Std. Error t-value   p-value    
#> b1:(Intercept)  -0.959953   0.098934 -9.7030 4.716e-11 ***
#> b2:(Intercept)  -1.414817   0.145057 -9.7535 4.160e-11 ***
#> d:(Intercept)  103.466985   2.772607 37.3176 < 2.2e-16 ***
#> e1:(Intercept)   9.209402   0.959596  9.5972 6.139e-11 ***
#> e1:(Intercept)   0.807378   0.072583 11.1236 1.574e-12 ***
#> f:(Intercept)    0.480612   0.273154  1.7595   0.08805 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  4.727843 (32 degrees of freedom)
```
