# bees

Data are from a binary mixture experiment that involves multiple
single-dose factorial designs where the insecticide imidacloprid is
combined with each of 7 pesticides in turn.

## Usage

``` r
data(bees)
```

## Format

A data frame with 66 observations on the following 5 variables.

- `mixture`:

  Indicator of single-dose experiment (or control)

- `treat`:

  Treatment or combination of treatments

- `rep`:

  Replication number (there were 3 replicates per group)

- `dead0h`:

  Number of dead bees initially

- `dead48h`:

  Number of dead bees after 48 hours

## Details

Imidacloprid is a widely used insecticide. In a recent study potential
synergistic effects on mortality of honey bees exposed to the insectide
in binary mixtures with seven pesticides from different classes:
acephate, λ-cyhalothrin, oxamyl, tetraconazole, sulfoxaflor, glyphosate,
and clothianidin were investigated. Bees were reared in cages (25
insects per cage), with three cages per treatment group, and exposed to
the mixture treatments for 48h. Mortality after 48h was the response.

## Source

Data were retrieved from PLoS ONE repository.

## References

Zhu YC, Yao J, Adamczyk J and Luttrell R, Synergistic toxicity and
physiological impact of imidacloprid alone and binary mixtures with
seven representative pesticides on honey bee (Apis mellifera). PLoS ONE
12: e0176837 (2017). https://doi.org/10.1371/journal.pone.0176837

## Examples

``` r
library(drc)

## Displaying the data
head(bees)
#>   mixture    treat rep dead0h dead48h
#> 1    CK-1      H2O   1      0       0
#> 2    CK-2      H2O   2      2       2
#> 3    CK-3      H2O   3      5       5
#> 4  AdvBra Advi274B   1      2       5
#> 5  AdvBra Advi274B   2      0       6
#> 6  AdvBra Advi274B   3      0       3

## Summarizing mortality by treatment
aggregate(dead48h ~ treat, data = bees, FUN = mean)
#>        treat   dead48h
#> 1   Adv274BL  4.666667
#> 2   Advi274B  4.666667
#> 3   Advi274D  4.666667
#> 4   Advi274K  4.666667
#> 5   Advi274R  4.666667
#> 6   Advi274T  4.666667
#> 7   Advi274V  4.666667
#> 8   AdviBela  8.000000
#> 9  AdviBrack 15.000000
#> 10  AdviDoma 10.666667
#> 11  AdviKara  5.666667
#> 12  AdviRoun  4.666667
#> 13  AdviTran  9.333333
#> 14  AdviVyda 15.000000
#> 15   Belay40  4.333333
#> 16   Brack91  9.333333
#> 17  Doma2500  2.666667
#> 18       H2O  2.333333
#> 19  Karat273  1.000000
#> 20  Roun2500  0.000000
#> 21  Trans117  2.000000
#> 22  Vydat162  6.666667
```
