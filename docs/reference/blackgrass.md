# Seedling Emergence of Blackgrass (Alopecurus myosuroides)

Seedling emergence of herbicide susceptible (S) and resistant (R)
Alopecurus myosuroides in reponse to sowing depth and suboptimal
temperature regimes (10/5C) and optimal temperature regimes (17/10C).

## Usage

``` r
data("blackgrass")
```

## Format

A data frame with 2752 observations on the following 12 variables.

- `Exp`:

  a numeric vector

- `Temp`:

  a numeric vector

- `Popu`:

  a numeric vector

- `Bio`:

  a factor with two levels

- `Depth`:

  a numeric vector

- `Rep`:

  a numeric vector

- `Start.Day`:

  a numeric vector

- `End.Day`:

  a numeric vector

- `Ger`:

  a numeric vector

- `Accum.Ger`:

  a numeric vector

- `TotalSeed`:

  a numeric vector

- `Pot`:

  a numeric vector

## References

Keshtkar, E., Mathiassen, S. K., Beffa, R., Kudsk, P. (2017). Seed
Germination and Seedling Emergence of Blackgrass (Alopecurus
myosuroides) as Affected by Non-Target-Site Herbicide Resistance. Weed
Science, 65, 732-742. https://doi.org/10.1017/wsc.2017.44

## Examples

``` r
library(drc)

## Displaying the data
head(blackgrass)
#>   Exp Temp Popu Bio Depth Rep Start.Day End.Day Ger Accum.Ger TotalSeed Pot
#> 1   1   10  914   S     0   1         0     360   0         0        36   3
#> 2   1   10  914   S     0   1       360     376   0         0        36   3
#> 3   1   10  914   S     0   1       376     384   0         0        36   3
#> 4   1   10  914   S     0   1       384     400   0         0        36   3
#> 5   1   10  914   S     0   1       400     408   0         0        36   3
#> 6   1   10  914   S     0   1       408     424   0         0        36   3

## Summarizing seedling emergence across treatments
aggregate(Accum.Ger ~ Temp + Bio, data = blackgrass, FUN = max)
#>   Temp Bio Accum.Ger
#> 1   10   R        33
#> 2   17   R        29
#> 3   10   S        36
#> 4   17   S        30
```
