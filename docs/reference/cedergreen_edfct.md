# Calculate Effective Dose for the Cedergreen-Ritz Hormesis Model

An internal helper function to calculate the effective dose (ED) and its
derivatives for the Cedergreen-Ritz five-parameter hormesis model. It
uses `uniroot` to find the dose for a given response level.

## Usage

``` r
cedergreen_edfct(
  parm,
  all_params,
  not_fixed,
  alpha,
  respl,
  reference,
  type,
  lower = 1e-04,
  upper = 10000
)
```

## Arguments

- parm:

  A numeric vector of the non-fixed model parameters.

- all_params:

  A numeric vector template for all model parameters (b,c,d,e,f).

- not_fixed:

  A logical or integer vector indicating the non-fixed parameters.

- alpha:

  A numeric value for the hormesis model's alpha shape parameter.

- respl:

  The response level to calculate the dose for (e.g., 50 for ED50).

- reference:

  A character string ("control" or "absolute") for calculating the
  response.

- type:

  A character string specifying the type of ED calculation.

- lower:

  The lower bound of the dose interval for the root-finding search.

- upper:

  The upper bound of the dose interval for the root-finding search.

## Value

A list containing the calculated effective dose and a vector of its
partial derivatives with respect to the non-fixed parameters.

## Author

Hannes Reinwald
