# Effect of technical grade and commercially formulated auxin herbicides

MCPA, 2,4-D, mecorprop and dichorlprop were applied either as technical
grades materials or as commercial formulations. Each experimental unit
consisted of five 1-week old seedlings grown together in a pot of
nutrient solution during 14 days.

## Usage

``` r
data(auxins)
```

## Format

A data frame with 150 observations on the following 5 variables.

- `dryweight`:

  a numeric vector

- `dose`:

  a numeric vector

- `replicate`:

  a factor with 3 levels

- `herbicide`:

  a factor with 5 levels

- `formulation`:

  a factor with 2 levels

## Details

Data are parts of a larger joint action experiment with various
herbicides.

The eight herbicide preparations are naturally grouped into four pairs
(herbicide:formulation) + control, and each pair of herbicides should
have the same active ingredients but different formulation constituents,
which were assumed to be biologically inert. The data consist of the 150
observations of dry weights, each observation being the weight of five
plants grown in the same pot. All the eight herbicide preparations have
essentially the same mode of action in the plant; they all act like the
plant auxins, which are plant regulators that affect cell enlongation an
other essential metabolic pathways. One of the objects of the experiment
was to test if the response functions were identical except for a
multiplicative factor in the dose. This is a necessary, but not a
sufficient, condition for a similar mode of action for the herbicides.

## Source

Streibig, J. C. (1987). Joint action of root-absorbed mixtures of auxin
herbicides in Sinapis alba L. and barley (Hordeum vulgare L.) *Weed
Research*, **27**, 337–347.

## References

Rudemo, M., Ruppert, D., and Streibig, J. C. (1989). Random-Effect
Models in Nonlinear Regression with Applications to Bioassay.
*Biometrics*, **45**, 349–362.

## Examples

``` r
library(drc)

## Displaying the data
head(auxins)
#>   dryweight  dose replicate herbicide formulation
#> 1      1.51 0.000         1   control     control
#> 2      1.43 0.000         1   control     control
#> 3      0.05 1.000         1      MCPA        tech
#> 4      0.06 0.500         1      MCPA        tech
#> 5      0.15 0.250         1      MCPA        tech
#> 6      0.40 0.125         1      MCPA        tech

## Fitting a four-parameter log-logistic model with different curves per herbicide
auxins.m1 <- drm(dryweight ~ dose, herbicide, data = auxins, fct = LL.4())
#> Control measurements detected for level: control
summary(auxins.m1)
#> 
#> Model fitted: Log-logistic (ED50 as parameter) (4 parms)
#> 
#> Parameter estimates:
#> 
#>                 Estimate Std. Error t-value   p-value    
#> b:MCPA         2.0707027  0.3893212  5.3188 4.261e-07 ***
#> b:24-D         2.5451091  0.7251643  3.5097 0.0006115 ***
#> b:mecorprop    2.1806515  0.6672591  3.2681 0.0013758 ** 
#> b:dichorlprop  1.3929604  0.6407011  2.1741 0.0314527 *  
#> c:MCPA         0.0476748  0.0481986  0.9891 0.3243816    
#> c:24-D         0.0267373  0.0545214  0.4904 0.6246531    
#> c:mecorprop    0.0714774  0.0669160  1.0682 0.2873660    
#> c:dichorlprop -0.0208598  0.3282497 -0.0635 0.9494243    
#> d:MCPA         1.2724185  0.0599127 21.2379 < 2.2e-16 ***
#> d:24-D         1.1472519  0.0877626 13.0722 < 2.2e-16 ***
#> d:mecorprop    1.2462043  0.1058507 11.7732 < 2.2e-16 ***
#> d:dichorlprop  1.2117451  0.1012722 11.9652 < 2.2e-16 ***
#> e:MCPA         0.0710631  0.0076089  9.3395 3.619e-16 ***
#> e:24-D         0.1275832  0.0147128  8.6716 1.241e-14 ***
#> e:mecorprop    0.1218995  0.0156854  7.7715 1.814e-12 ***
#> e:dichorlprop  0.3391793  0.1392770  2.4353 0.0161949 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  0.1535001 (134 degrees of freedom)

## Plotting the fitted curves
plot(auxins.m1, xlab = "Dose", ylab = "Dry weight")
```
