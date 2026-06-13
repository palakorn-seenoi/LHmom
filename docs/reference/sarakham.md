# Annual Maximum Daily Temperature for Maha Sarakham, Thailand

A dataset containing the annual maximum series (AMS) of daily
temperature measurements in Maha Sarakham, Thailand, from 1985 to 2025.

## Usage

``` r
data(sarakham)
```

## Format

A data frame with 41 rows and 2 variables:

- year:

  The observation year (1985-2025).

- temperature:

  The annual maximum daily temperature in degrees Celsius (°C).

## Source

Thai Meteorological Department (TMD)

## Examples

``` r
data(sarakham)
head(sarakham)
#>   year temperature
#> 1 1985        41.2
#> 2 1986        41.0
#> 3 1987        41.0
#> 4 1988        40.5
#> 5 1989        39.5
#> 6 1990        40.3
```
