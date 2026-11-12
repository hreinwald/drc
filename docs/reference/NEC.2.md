# Two-parameter NEC model

Convenience function for the NEC model with lower limit fixed at 0 and
upper limit fixed.

## Usage

``` r
NEC.2(upper = 1, fixed = c(NA, NA), names = c("b", "e"), ...)
```

## Arguments

- upper:

  numeric value. The fixed upper limit in the model. Default is 1.

- fixed:

  numeric vector of length 2 specifying fixed parameters (NAs for free
  parameters).

- names:

  character vector of parameter names.

- ...:

  additional arguments passed to
  [`NEC`](https://hreinwald.github.io/drc/reference/NEC.md).

## Value

A list (see [`NEC`](https://hreinwald.github.io/drc/reference/NEC.md)).

## See also

[`NEC`](https://hreinwald.github.io/drc/reference/NEC.md),
[`NEC.3`](https://hreinwald.github.io/drc/reference/NEC.3.md),
[`NEC.4`](https://hreinwald.github.io/drc/reference/NEC.4.md)

## Examples

``` r
## Response relative to the control mean, upper limit fixed at 1
ryegrass.rel <- transform(ryegrass, rel = rootl / mean(rootl[conc == 0]))
nec.m2 <- drm(rel ~ conc, data = ryegrass.rel, fct = NEC.2(upper = 1))
plot(nec.m2, type = "all", main = "NEC.2 (c = 0, d = 1)")

```
