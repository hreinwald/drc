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
#>           1          2         3          4          5          6
#> 5 1.8123778 0.70797936 0.7069777 0.74315403 0.03049114 0.47266280
#> 6 0.7690992 0.07492129 0.1991684 0.50155244 0.02740665 0.26126638
#> 7 0.6766928 0.53788677 0.3238083 0.57279347 0.49562667 0.02968098
#> 8 0.5990200 0.20800038 0.2113740 0.05549486 0.25613766 0.28737209
#> 9 0.1800618 0.76738087 0.5514389 0.11201090 0.48343071 0.46126469
#> 
#> 
#> ED value considered: 50 
#> Conc. no.\Replicates: 
#>           1          2          3          4         5         6
#> 5 0.6991675 1.37235074 0.49271784 0.44941966 0.1074398 0.6018120
#> 6 7.3839802 0.39097648 0.09645955 0.58216624 0.4309672 0.4252722
#> 7 0.9153291 0.02106256 0.46655957 0.20210381 0.3129893 0.3744315
#> 8 0.3726099 0.16439252 0.21873962 0.08030052 0.2577437 0.2953281
#> 9 0.3519223 1.19099163 0.34465510 0.13020887 0.4812851 0.1996946
#> 
#> 
```
