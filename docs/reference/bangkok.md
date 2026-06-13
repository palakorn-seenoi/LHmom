# Annual Maximum Daily Rainfall for Bangkok, Thailand

A dataset containing the annual maximum series (AMS) of daily rainfall
measurements in Bangkok, Thailand, from 1951 to 2025. This dataset is
typically used for extreme value analysis and hydrological frequency
modeling.

## Usage

``` r
data(bangkok)
```

## Format

A data frame with 75 rows and 2 variables:

- year:

  The observation year (1951-2025).

- rainfall:

  The annual maximum daily rainfall in millimeters (mm).

## Source

Thai Meteorological Department (TMD)

## Examples

``` r
data(bangkok)
head(bangkok)
#>   year rainfall
#> 1 1951    133.5
#> 2 1952    111.0
#> 3 1953     84.1
#> 4 1954     53.8
#> 5 1955    108.8
#> 6 1956     69.4
```
