# Extracting residuals from the fitted dose-response model

`residuals` extracts different types of residuals from an object of
class 'drc'.

## Usage

``` r
# S3 method for class 'drc'
residuals(
  object,
  typeRes = c("working", "standardised", "studentised"),
  trScale = TRUE,
  ...
)
```

## Arguments

- object:

  an object of class 'drc'.

- typeRes:

  character string specifying the type of residual to be returned:
  raw/working residuals, residuals standardised using the estimated
  residual standard error, or studentised residuals based on the H
  matrix of partial derivatives of the model function.

- trScale:

  logical value indicating whether or not to return residuals on the
  transformed scale (in case a Box-Cox transformation was applied).

- ...:

  additional arguments.

## Value

The raw (also called working) residuals or some kind of scaled residuals
extracted from `object`.

## Author

Christian Ritz

## Examples

``` r
## Fitting a four-parameter log-logistic model
ryegrass.m1 <- drm(rootl ~ conc, data = ryegrass, fct = LL.4())

## Displaying the residual plot (raw residuals)
plot(fitted(ryegrass.m1), residuals(ryegrass.m1))


## Using the standardised residuals
plot(fitted(ryegrass.m1), residuals(ryegrass.m1, typeRes = "standard"))

```
