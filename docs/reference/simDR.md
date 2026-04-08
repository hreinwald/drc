# Simulating ED values under various scenarios

Simulating ED values for a given model and given dose values.

## Usage

``` r
simDR(
  mpar,
  sigma,
  fct,
  noSim = 1000,
  conc,
  edVec = c(10, 50),
  seedVal = 20070723
)
```

## Arguments

- mpar:

  numeric vector of model parameters.

- sigma:

  numeric specifying the residual standard deviation.

- fct:

  list supplying the chosen dose-response mean function (e.g.,
  [`LL.4()`](https://hreinwald.github.io/drc/reference/LL.4.md)).

- noSim:

  numeric giving the number of simulations. Defaults to `1000`.

- conc:

  numeric vector of concentration/dose values. Must contain at least 5
  values.

- edVec:

  numeric vector of ED levels to estimate in each simulation. Defaults
  to `c(10, 50)`.

- seedVal:

  numeric giving the seed used to initialise the random number
  generator. Defaults to `20070723`.

## Value

Invisibly returns a list with one element:

- `se`:

  A 3D array of dimensions `(length(conc) - 4) x 6 x length(edVec)`
  containing empirical standard deviations of the estimated ED values.
  Rows correspond to the number of concentration levels used (starting
  from 5). Columns correspond to the number of replicates per
  concentration (1 to 6). The third dimension corresponds to each ED
  level in `edVec`.

The array values are also printed to the console during execution.

## Details

The arguments `mpar` and `sigma` are typically obtained from a previous
model fit. Only dose-response models assuming normally distributed
errors can be used.

## Author

Christian Ritz, Hannes Reinwald

## Examples

``` r
ryegrass.m1 <- drm(ryegrass, fct = LL.4())

simDR(
  mpar     = coef(ryegrass.m1),
  sigma    = sqrt(summary(ryegrass.m1)$resVar),
  fct      = LL.4(),
  noSim    = 2,
  conc     = c(1.88, 3.75, 7.50, 0.94, 15, 0.47, 30, 0.23, 60),
  seedVal  = 20070723
)
#> Concentrations used: 1.88 3.75 7.5 0.94 15 0.47 30 0.23 60 
#> 
#> ED value considered: 10 
#> Conc. no.\Replicates: 
#>           1          2         3          4         5          6
#> 5 1.7937943 0.70676461 0.7070331 0.74477323 0.0314240 0.47265995
#> 6 0.7715370 0.07489298 0.1991687 0.50345455 0.0274019 0.26133303
#> 7 0.6766681 0.43816932 0.3241948 0.57191204 0.4956274 0.02968537
#> 8 0.5276024 0.20800084 0.2113615 0.05549555 0.2561400 0.28737212
#> 9 0.1800602 0.76739811 0.5514185 0.11201042 0.4834278 0.46828149
#> 
#> 
#> ED value considered: 50 
#> Conc. no.\Replicates: 
#>           1          2          3         4         5         6
#> 5 0.6998346 1.32865730 0.49266733 0.4490645 0.1069289 0.6018208
#> 6 6.7983515 0.39097416 0.09645987 0.5823232 0.4309704 0.4252761
#> 7 0.9155897 0.03037852 0.46647071 0.2017280 0.3129933 0.3744331
#> 8 0.3950825 0.16439311 0.21879117 0.0803034 0.2577441 0.2953281
#> 9 0.3519176 1.19098825 0.34427430 0.1302080 0.4812865 0.2009490
#> 
#> 
```
