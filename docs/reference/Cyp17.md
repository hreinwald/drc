# Cyp17 expression data

Observed Cyp17 gene expression measured at several dose levels across
multiple experimental runs. CYP17 is a key enzyme in steroid hormone
biosynthesis, and changes in its expression can indicate
endocrine-disrupting effects.

## Usage

``` r
data(Cyp17)
```

## Format

A data frame with 63 observations on the following 3 variables.

- `run`:

  ID of 3 different runs

- `dose`:

  5 dose levels (0, 0.1, 10, 100, 500)

- `expression`:

  observed expression

## Examples

``` r
data(Cyp17)

## Display the structure of the data
head(Cyp17)
#>   run dose expression
#> 1   1  500 0.04934222
#> 2   1  500 0.05863851
#> 3   1  500 0.07518652
#> 4   1  100 0.05492933
#> 5   1  100 0.03329362
#> 6   1  100 0.01662889

## Log-transform the expression values
Cyp17$logexpression <- log(Cyp17$expression) + 5

## Fit a four-parameter log-logistic model (ignoring run effects)
Cyp17.m1 <- drm(logexpression ~ dose, data = Cyp17, fct = LL.4())
summary(Cyp17.m1)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) (4 parms)
#> 
#> Parameter estimates:
#> 
#>                 Estimate Std. Error t-value   p-value    
#> b:(Intercept)   -0.53357    0.49317 -1.0819    0.2837    
#> c:(Intercept)    1.07435    0.10463 10.2682 9.388e-15 ***
#> d:(Intercept)    2.49393    2.12401  1.1742    0.2450    
#> e:(Intercept)  556.15323 2932.24810  0.1897    0.8502    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.5298351 (59 degrees of freedom)
plot(Cyp17.m1, main = "Cyp17 dose-response")
```
