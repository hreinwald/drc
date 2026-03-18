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

Christian Ritz

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
#> Michaelis-Menten 
#> (2 parameters) 
#> In 'drc':  MM.2 
#> 
#> Shifted Michaelis-Menten 
#> (3 parameters) 
#> In 'drc':  MM.3 
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

## Listing all (log-)logistic functions
getMeanFunctions(fname = "L")

## Listing all three-parameter (log-)logistic or Weibull functions
getMeanFunctions(3, fname = c("LL", "W"))
```
