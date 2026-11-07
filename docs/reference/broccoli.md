# The Effects of Drought Stress on Leaf Development in a *Brassica oleracea* population

The effect of drought stress on *Brassica oleracea* should be
investigated, selecting drought stress resistant out of a population of
different DH genotypes. The study was carried out on 48 DH lines
developed from F1 plants of a cross between the rapid cycling chinese
kale (*Brassica oleracea* var. *alboglabra* (L.H. Bailey) Musil) and
broccoli (*Brassica oleracea* var. *italica* Plenck). 2 stress
treatments (not watered and a watered control) are randomly assigned to
4 plants per genotype (2 per treatment) resulting in 192 plants in
total. For the genotypes 5, 17, 31, 48, additional 12 plants (6 per
treatment) are included into the completely randomized design, which
results in a total of 240 plants. For each plant the length of the
youngest leaf at the beginning of the experiment is measured daily for a
period of 16 days. For the additional 12 plants of the 4 genotypes the
leaf water potential was measured as a secondary endpoint (omitted
here); due to these destructive measurements some dropouts occur.

## Usage

``` r
data(broccoli)
```

## Format

A data frame with 3689 observations on the following 5 variables.

- `LeafLength`:

  Length of the youngest leaf \[cm\]

- `ID`:

  Plant identifier for 240 plants

- `Stress`:

  Drought stress treatment with 2 levels (control/drought)

- `Genotype`:

  Genotype ID with 48 levels

- `Day`:

  Day of repeated measurement (1,2,...,16)

## References

Uptmoor, R.; Osei-Kwarteng, M.; Guertler, S. & Stuetzel, H. Modeling the
Effects of Drought Stress on Leaf Development in a Brassica oleracea
Doubled Haploid Population Using Two-phase Linear Functions. Journal of
the American Society for Horticultural Science, 2009, 134, 543-552.

## Examples

``` r
data(broccoli)

## Display the structure of the data
head(broccoli)
#>   LeafLength  ID  Stress Genotype Day
#> 1        1.4  38 control       17   1
#> 2        1.2  62 control       17   1
#> 3        2.5  35 control       17   1
#> 4        1.8  91 control       17   1
#> 5        1.7  76 control       17   1
#> 6        1.4 108 control       17   1

## Fit a five-parameter log-logistic model per stress treatment
broccoli.m1 <- drm(LeafLength ~ Day, curveid = Stress,
                   data = broccoli, fct = LL.5())
summary(broccoli.m1)
#> 
#> Model fitted: Generalized log-logistic (ED50 as parameter) (5 parms)
#> 
#> Parameter estimates:
#> 
#>           Estimate Std. Error t-value   p-value    
#> b:control -3.81852    0.81127 -4.7068 2.609e-06 ***
#> b:drought -4.63724    1.03424 -4.4837 7.560e-06 ***
#> c:control  1.87071    0.30308  6.1723 7.463e-10 ***
#> c:drought  2.02890    0.38197  5.3117 1.151e-07 ***
#> d:control 14.58605    0.38050 38.3342 < 2.2e-16 ***
#> d:drought 10.88021    0.14644 74.2980 < 2.2e-16 ***
#> e:control  9.62040    0.82948 11.5981 < 2.2e-16 ***
#> e:drought  7.97475    0.72091 11.0620 < 2.2e-16 ***
#> f:control  0.38222    0.13522  2.8267  0.004729 ** 
#> f:drought  0.28037    0.11034  2.5410  0.011094 *  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Residual standard error:
#> 
#>  1.69312 (3679 degrees of freedom)
plot(broccoli.m1, main = "Broccoli leaf growth by stress treatment")
```
