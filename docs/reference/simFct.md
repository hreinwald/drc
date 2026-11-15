# Simulation of dose-response data and ED estimation

Simulates dose-response datasets using parametric or non-parametric
methods and estimates effective doses (ED values) from each simulated
dataset. Useful for assessing the performance of ED estimation methods
via Monte Carlo simulation.

## Usage

``` r
simFct(
  noSim,
  edVal = c(10, 20, 50),
  type = c("non-parametric", "parametric"),
  response = c("bin", "con"),
  fct = LL.2(),
  coefVec,
  method = c("sp", "p", "np"),
  doseVec,
  nVec,
  pVec,
  rVec,
  resVar,
  pfct = fct,
  reference = NULL,
  span = NA,
  minmax = "response",
  lower = NULL,
  upper = NULL,
  seedVal = 200810201
)
```

## Arguments

- noSim:

  integer. Number of simulations to run.

- edVal:

  numeric vector of ED levels to estimate (default is `c(10, 20, 50)`).

- type:

  character string. Either "non-parametric" or "parametric" simulation.

- response:

  character string. Either "bin" (binomial) or "con" (continuous)
  response.

- fct:

  dose-response function used for simulation (default is
  [`LL.2()`](https://hreinwald.github.io/drc/reference/LL.2.md)).

- coefVec:

  numeric vector of model coefficients for parametric simulation.

- method:

  character string. Estimation method: "sp" (semi-parametric), "p"
  (parametric), or "np" (non-parametric).

- doseVec:

  numeric vector of dose values.

- nVec:

  numeric vector of sample sizes per dose (for binomial response).

- pVec:

  numeric vector of expected response probabilities (for non-parametric
  simulation).

- rVec:

  numeric vector of responses.

- resVar:

  numeric. Residual variance (for continuous response).

- pfct:

  dose-response function used for fitting (defaults to `fct`).

- reference:

  character string specifying the reference for ED estimation.

- span:

  numeric. Smoothing parameter for local regression. NA uses default.

- minmax:

  character string. Type of min/max calculation. Default is "response".

- lower:

  numeric. Lower bounds for optimization.

- upper:

  numeric. Upper bounds for optimization.

- seedVal:

  integer. Random seed for reproducibility (default is 200810201).

## Value

A list with components `edArray` (array of ED estimates), `mixVec`,
`edVal`, `aicVec`, and `spanVec`.

## Author

Christian Ritz
