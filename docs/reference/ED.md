# Estimating effective doses

S3 generic function that dispatches to the appropriate method for
estimating effective concentrations (EC) or effective doses (ED) at
specified response levels. For objects of class `drc`, the default
method [`ED.drc`](https://hreinwald.github.io/drc/reference/ED.drc.md)
is called.

## Usage

``` r
ED(object, ...)
```

## Arguments

- object:

  an object of class `drc`.

- ...:

  additional arguments passed to the method.

## Value

See [`ED.drc`](https://hreinwald.github.io/drc/reference/ED.drc.md) for
details on the return value.

## See also

[`ED.drc`](https://hreinwald.github.io/drc/reference/ED.drc.md) for the
default method,
[`EDcomp`](https://hreinwald.github.io/drc/reference/EDcomp.md) for
estimating differences and ratios of ED

## Author

Christian Ritz
