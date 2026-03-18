# Convert absolute to relative response levels

Internal helper that converts an absolute response level to a relative
(percentage) scale based on the upper and lower asymptotes of a
dose-response curve.

## Usage

``` r
absToRel(parmVec, respl, typeCalc)
```

## Arguments

- parmVec:

  numeric vector of model parameters where the third element is the
  upper asymptote and the second element is the lower asymptote.

- respl:

  numeric response level to convert.

- typeCalc:

  character string. If "absolute", the conversion is performed;
  otherwise the input `respl` is returned unchanged.

## Value

A numeric value representing the (possibly converted) response level as
a percentage.
