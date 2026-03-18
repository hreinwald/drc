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

  list supplying the chosen mean function.

- noSim:

  numeric giving the number of simulations.

- conc:

  numeric vector of concentration/dose values.

- edVec:

  numeric vector of ED values to estimate in each simulation.

- seedVal:

  numeric giving the seed used to initiate the random number generator.

## Value

A list of matrices with as many components as there are chosen ED
values. The entries in the matrices are empirical standard deviations of
the estimated ED values. Row-wise from top to bottom more and more
concentration/dose values are included in the simulations; top row
starting with 5 concentrations. The number of replicates increases
column by column from left to right.

The list is returned invisibly as the matrices also are displayed.

## Details

The arguments `mpar` and `sigma` are typically obtained from a previous
model fit.

Only dose-response models assuming normally distributed errors can be
used.

## Author

Christian Ritz

## Examples

``` r
ryegrass.m1 <- drm(ryegrass, fct=LL.4())

simDR(coef(ryegrass.m1), sqrt(summary(ryegrass.m1)$resVar), LL.4(), 2,
c(1.88, 3.75, 7.50, 0.94, 15, 0.47, 30, 0.23, 60), seedVal = 200710291)
#> Concentrations used: 1.88 3.75 7.5 0.94 15 0.47 30 0.23 60 
#> 
#> ED value considered: 10 
#> Conc. no.\Replicates: 
#>            1          2          3          4          5          6
#> 5 0.21663344 0.81050689 0.03720474 0.31941962 0.27602174 1.06748994
#> 6 0.94123238 0.09153931 0.05169959 0.42236519 0.05903568 0.65652040
#> 7 0.37505197 0.89749084 1.10562748 0.46717540 0.09427895 0.09941082
#> 8 0.29322754 0.35947683 0.36292195 0.12434167 0.25382652 0.44720085
#> 9 0.05051407 0.81151253 0.36009230 0.09176809 0.05572999 0.07993389
#> 
#> 
#> ED value considered: 50 
#> Conc. no.\Replicates: 
#>            1         2         3         4         5         6
#> 5 0.41931241 2.5598699 0.1999825 0.1641816 0.0184324 1.2226656
#> 6 0.59150811 0.1644054 0.1555681 0.3922277 0.1010897 0.6564398
#> 7 0.79504087 0.5011138 0.4209378 0.3844431 0.3784710 0.4346236
#> 8 0.01548633 0.4687979 0.4580110 0.6642111 0.2163613 0.1593073
#> 9 1.82422377 0.4866429 0.5804128 0.1243123 0.2344158 0.1205998
#> 
#> 
```
