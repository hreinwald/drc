# Display available dose-response models

Display information about available, built-in dose-response models. The
arguments `noParm` and `fname` can be combined.

## Usage

``` r
getMeanFunctions(noParm = NA, fname = NULL, flist = NULL, display = TRUE)
```

## Arguments

- noParm:

  numeric specifying the number of parameters of the models to be
  displayed. The default (NA) results in display of all models,
  regardless of number of parameters.

- fname:

  character string or vector of character strings specifying the short
  name(s) of the models to be displayed (need to match exactly).

- flist:

  list of built-in functions to be displayed.

- display:

  logical indicating whether or not the requested models should be
  displayed on the R console.

## Value

An invisible list of functions or a list of strings with brief function
descriptions.

## Author

Christian Ritz, Hannes Reinwald

## Examples

``` r
## Listing all functions
getMeanFunctions()
#> Log-logistic (ED50 as parameter) with lower limit at 0 and upper limit at 1 
#> (2 parameters) 
#> In 'drc':  LL.2 
#> 
#> Log-logistic (ED50 as parameter) with lower limit at 0 
#> (3 parameters) 
#> In 'drc':  LL.3 
#> 
#> Log-logistic (ED50 as parameter) with upper limit at 1 
#> (3 parameters) 
#> In 'drc':  LL.3u 
#> 
#> Log-logistic (ED50 as parameter) 
#> (4 parameters) 
#> In 'drc':  LL.4 
#> 
#> Generalized log-logistic (ED50 as parameter) 
#> (5 parameters) 
#> In 'drc':  LL.5 
#> 
#> Weibull (type 1) with lower limit at 0 and upper limit at 1 
#> (2 parameters) 
#> In 'drc':  W1.2 
#> 
#> Weibull (type 1) with lower limit at 0 
#> (3 parameters) 
#> In 'drc':  W1.3 
#> 
#> Weibull (type 1) with upper limit at 1 
#> (3 parameters) 
#> In 'drc':  W1.3u 
#> 
#> Weibull (type 1) 
#> (4 parameters) 
#> In 'drc':  W1.4 
#> 
#> Weibull (type 2) with lower limit at 0 and upper limit at 1 
#> (2 parameters) 
#> In 'drc':  W2.2 
#> 
#> Weibull (type 2) with lower limit at 0 
#> (3 parameters) 
#> In 'drc':  W2.3 
#> 
#> Weibull (type 2) with upper limit at 1 
#> (3 parameters) 
#> In 'drc':  W2.3u 
#> 
#> Weibull (type 2) 
#> (4 parameters) 
#> In 'drc':  W2.4 
#> 
#> Brain-Cousens (hormesis) with lower limit fixed at 0 
#> (4 parameters) 
#> In 'drc':  BC.4 
#> 
#> Brain-Cousens (hormesis) 
#> (5 parameters) 
#> In 'drc':  BC.5 
#> 
#> Log-logistic (log(ED50) as parameter) with lower limit at 0 and upper limit at 1 
#> (2 parameters) 
#> In 'drc':  LL2.2 
#> 
#> Log-logistic (log(ED50) as parameter) with lower limit at 0 
#> (3 parameters) 
#> In 'drc':  LL2.3 
#> 
#> Log-logistic (log(ED50) as parameter) with upper limit at 1 
#> (3 parameters) 
#> In 'drc':  LL2.3u 
#> 
#> Log-logistic (log(ED50) as parameter) 
#> (4 parameters) 
#> In 'drc':  LL2.4 
#> 
#> Generalised log-logistic (log(ED50) as parameter) 
#> (5 parameters) 
#> In 'drc':  LL2.5 
#> 
#> Asymptotic regression with lower limit at 0 
#> (2 parameters) 
#> In 'drc':  AR.2 
#> 
#> Shifted asymptotic regression 
#> (3 parameters) 
#> In 'drc':  AR.3 
#> 
#> Exponential decay with lower limit at 0 
#> (2 parameters) 
#> In 'drc':  EXD.2 
#> 
#> Shifted exponential decay 
#> (3 parameters) 
#> In 'drc':  EXD.3 
#> 
#> Gompertz with lower limit at 0 and upper limit at 1 
#> (2 parameters) 
#> In 'drc':  G.2 
#> 
#> Gompertz with lower limit at 0 
#> (3 parameters) 
#> In 'drc':  G.3 
#> 
#> Gompertz with upper limit at 1 
#> (3 parameters) 
#> In 'drc':  G.3u 
#> 
#> Gompertz 
#> (4 parameters) 
#> In 'drc':  G.4 
#> 
#> Log-normal with lower limit at 0 and upper limit at 1 
#> (2 parameters) 
#> In 'drc':  LN.2 
#> 
#> Log-normal with lower limit at 0 
#> (3 parameters) 
#> In 'drc':  LN.3 
#> 
#> Log-normal with upper limit at 1 
#> (3 parameters) 
#> In 'drc':  LN.3u 
#> 
#> Log-normal 
#> (4 parameters) 
#> In 'drc':  LN.4 
#> 
#> Logistic (ED50 as parameter) with lower limit fixed at 0 
#> (3 parameters) 
#> In 'drc':  L.3 
#> 
#> Logistic (ED50 as parameter) 
#> (4 parameters) 
#> In 'drc':  L.4 
#> 
#> Generalised logistic (ED50 as parameter) 
#> (5 parameters) 
#> In 'drc':  L.5 
#> 
#> Michaelis-Menten 
#> (2 parameters) 
#> In 'drc':  MM.2 
#> 
#> Shifted Michaelis-Menten 
#> (3 parameters) 
#> In 'drc':  MM.3 
#> 
#> Cedergreen-Ritz-Streibig with lower limit 0 (alpha=1) 
#> (4 parameters) 
#> In 'drc':  CRS.4a 
#> 
#> Cedergreen-Ritz-Streibig with lower limit 0 (alpha=0.5) 
#> (4 parameters) 
#> In 'drc':  CRS.4b 
#> 
#> Cedergreen-Ritz-Streibig with lower limit 0 (alpha=0.25) 
#> (4 parameters) 
#> In 'drc':  CRS.4c 
#> 
#> Cedergreen-Ritz-Streibig (alpha=1) 
#> (5 parameters) 
#> In 'drc':  CRS.5a 
#> 
#> Cedergreen-Ritz-Streibig (alpha=0.5) 
#> (5 parameters) 
#> In 'drc':  CRS.5b 
#> 
#> Cedergreen-Ritz-Streibig (alpha=0.25) 
#> (5 parameters) 
#> In 'drc':  CRS.5c 
#> 
#> Generalised Cedergreen-Ritz-Streibig (hormesis) 
#> (6 parameters) 
#> In 'drc':  CRS.6 
#> 
#> U-shaped Cedergreen-Ritz-Streibig 
#> (4 parameters) 
#> In 'drc':  UCRS.4a 
#> 
#> U-shaped Cedergreen-Ritz-Streibig 
#> (4 parameters) 
#> In 'drc':  UCRS.4b 
#> 
#> U-shaped Cedergreen-Ritz-Streibig 
#> (4 parameters) 
#> In 'drc':  UCRS.4c 
#> 
#> U-shaped Cedergreen-Ritz-Streibig 
#> (5 parameters) 
#> In 'drc':  UCRS.5a 
#> 
#> U-shaped Cedergreen-Ritz-Streibig 
#> (5 parameters) 
#> In 'drc':  UCRS.5b 
#> 
#> U-shaped Cedergreen-Ritz-Streibig 
#> (5 parameters) 
#> In 'drc':  UCRS.5c 
#> 
#> NEC with lower limit at 0 and upper limit at 1 
#> (2 parameters) 
#> In 'drc':  NEC.2 
#> 
#> NEC with lower limit at 0 
#> (3 parameters) 
#> In 'drc':  NEC.3 
#> 
#> NEC 
#> (4 parameters) 
#> In 'drc':  NEC.4 
#> 
#> Fractional polynomial 
#> (4 parameters) 
#> In 'drc':  FPL.4(-1,1) 
#> 
#> Weibull (type 2) with lower limit at 0 
#> (3 parameters) 
#> In 'drc':  W2x.3 
#> 
#> Weibull (type 2) 
#> (4 parameters) 
#> In 'drc':  W2x.4 
#> 
#> Baroreflex 
#> (5 parameters) 
#> In 'drc':  baro5 
#> 
#> Gamma 
#> (4 parameters) 
#> In 'drc':  gammadr 
#> 

## Listing all functions with 4 parameters
getMeanFunctions(4)
#> Log-logistic (ED50 as parameter) 
#> (4 parameters) 
#> In 'drc':  LL.4 
#> 
#> Weibull (type 1) 
#> (4 parameters) 
#> In 'drc':  W1.4 
#> 
#> Weibull (type 2) 
#> (4 parameters) 
#> In 'drc':  W2.4 
#> 
#> Brain-Cousens (hormesis) with lower limit fixed at 0 
#> (4 parameters) 
#> In 'drc':  BC.4 
#> 
#> Log-logistic (log(ED50) as parameter) 
#> (4 parameters) 
#> In 'drc':  LL2.4 
#> 
#> Gompertz 
#> (4 parameters) 
#> In 'drc':  G.4 
#> 
#> Log-normal 
#> (4 parameters) 
#> In 'drc':  LN.4 
#> 
#> Logistic (ED50 as parameter) 
#> (4 parameters) 
#> In 'drc':  L.4 
#> 
#> Cedergreen-Ritz-Streibig with lower limit 0 (alpha=1) 
#> (4 parameters) 
#> In 'drc':  CRS.4a 
#> 
#> Cedergreen-Ritz-Streibig with lower limit 0 (alpha=0.5) 
#> (4 parameters) 
#> In 'drc':  CRS.4b 
#> 
#> Cedergreen-Ritz-Streibig with lower limit 0 (alpha=0.25) 
#> (4 parameters) 
#> In 'drc':  CRS.4c 
#> 
#> U-shaped Cedergreen-Ritz-Streibig 
#> (4 parameters) 
#> In 'drc':  UCRS.4a 
#> 
#> U-shaped Cedergreen-Ritz-Streibig 
#> (4 parameters) 
#> In 'drc':  UCRS.4b 
#> 
#> U-shaped Cedergreen-Ritz-Streibig 
#> (4 parameters) 
#> In 'drc':  UCRS.4c 
#> 
#> NEC 
#> (4 parameters) 
#> In 'drc':  NEC.4 
#> 
#> Fractional polynomial 
#> (4 parameters) 
#> In 'drc':  FPL.4(-1,1) 
#> 
#> Weibull (type 2) 
#> (4 parameters) 
#> In 'drc':  W2x.4 
#> 
#> Gamma 
#> (4 parameters) 
#> In 'drc':  gammadr 
#> 

## Listing all (log-)logistic functions
getMeanFunctions(fname = "L")

## Listing all three-parameter (log-)logistic or Weibull functions
getMeanFunctions(3, fname = c("LL", "W"))
```
