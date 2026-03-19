# Extract fitted values from model

Extracts fitted values from an object of class 'drc'.

## Usage

``` r
# S3 method for class 'drc'
fitted(object, ...)
```

## Arguments

- object:

  an object of class 'drc'.

- ...:

  additional arguments.

## Value

Fitted values extracted from `object`.

## Author

Christian Ritz

## Examples

``` r
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())
plot(fitted(ryegrass.m1), residuals(ryegrass.m1))  # a residual plot

```
