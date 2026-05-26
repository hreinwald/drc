# Sets control arguments

Set control arguments in the control argument in the function
[`drm`](https://hreinwald.github.io/drc/reference/drm.md).

## Usage

``` r
drmc(
  constr = FALSE,
  errorm = TRUE,
  maxIt = 500,
  method = "BFGS",
  noMessage = FALSE,
  relTol = 1e-10,
  rmNA = FALSE,
  useD = FALSE,
  trace = FALSE,
  otrace = FALSE,
  warnVal = -1,
  dscaleThres = 1e-15,
  rscaleThres = 1e-15,
  conCheck = TRUE
)
```

## Arguments

- constr:

  logical. If `TRUE` optimisation is constrained, only yielding
  non-negative parameters.

- errorm:

  logical specifying whether failed convergence in
  [`drm`](https://hreinwald.github.io/drc/reference/drm.md) should
  result in an error or only a warning.

- maxIt:

  numeric. The maximum number of iterations in the optimisation
  procedure.

- method:

  character string. The method used in the optimisation procedure. See
  [`optim`](https://rdrr.io/r/stats/optim.html) for available methods.

- noMessage:

  logical, specifying whether or not messages should be displayed.

- relTol:

  numeric. The relative tolerance in the optimisation procedure. A
  tighter tolerance (smaller value) improves cross-platform
  reproducibility of results by ensuring the optimiser converges closer
  to the true optimum regardless of platform-specific floating-point
  behaviour. Default is `1e-10`.

- rmNA:

  logical. Should `NA`s be removed from sum of squares used for
  estimation? Default is `FALSE` (not removed).

- useD:

  logical. If `TRUE` derivatives are used for estimation (if available).

- trace:

  logical. If `TRUE` the trace from
  [`optim`](https://rdrr.io/r/stats/optim.html) is displayed.

- otrace:

  logical. If `TRUE` error messages from the optimisation are displayed.

- warnVal:

  numeric. If equal to 0 then the warnings are stored and displayed at
  the end. See under ‘warn’ in
  [`options`](https://rdrr.io/r/base/options.html). The default results
  in suppression of warnings.

- dscaleThres:

  numeric value specifying the threshold for dose scaling.

- rscaleThres:

  numeric value specifying the threshold for response scaling.

- conCheck:

  logical, switching on/off handling of control measurements.

## Value

A list with components corresponding to each of the above arguments.

## See also

[`drm`](https://hreinwald.github.io/drc/reference/drm.md),
[`optim`](https://rdrr.io/r/stats/optim.html)

## Author

Christian Ritz

## Examples

``` r
## Displaying the default settings
drmc()
#> $constr
#> [1] FALSE
#> 
#> $errorm
#> [1] TRUE
#> 
#> $maxIt
#> [1] 500
#> 
#> $method
#> [1] "BFGS"
#> 
#> $noMessage
#> [1] FALSE
#> 
#> $relTol
#> [1] 1e-07
#> 
#> $rmNA
#> [1] FALSE
#> 
#> $useD
#> [1] FALSE
#> 
#> $trace
#> [1] FALSE
#> 
#> $otrace
#> [1] FALSE
#> 
#> $warnVal
#> [1] -1
#> 
#> $dscaleThres
#> [1] 1e-15
#> 
#> $rscaleThres
#> [1] 1e-15
#> 
#> $conCheck
#> [1] TRUE
#> 

## Using the 'method' argument
model1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
model2 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4(),
  control = drmc(method = "Nelder-Mead"))
```
