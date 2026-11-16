# Effective dose for the U-shaped Cedergreen-Ritz-Streibig model

Internal helper returning the ED and its gradient for
[`ucedergreen`](https://hreinwald.github.io/drc/reference/ucedergreen.md).
The U-shaped curve starts at `c`, dips, and then increases towards `d`
(`b > 0`). The dose is searched on the increasing branch beyond the dip,
where the equation has a unique root.

## Usage

``` r
ucedergreen_edfct(
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

  numeric vector of the non-fixed parameters.

- all_params:

  numeric vector template with all parameters (b, c, d, e, f).

- not_fixed:

  logical vector flagging the non-fixed parameters.

- alpha:

  the fixed shape parameter.

- respl:

  response level: percent of the way from `c` to `d` for
  `type = "relative"`, the response itself for `type = "absolute"`.

- reference:

  `"control"` (percentage measured from `c`) or otherwise percentage
  measured down from `d`.

- type:

  `"relative"` or `"absolute"`.

- lower, upper:

  dose search interval.

## Value

list(ED, gradient of ED w.r.t. the non-fixed parameters); `NA`s with a
warning if the level is not reached.
