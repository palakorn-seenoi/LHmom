# Calculate the Lambda Function for the PE3 Distribution

An internal helper function to compute the weighted sum of the Omega
function values, calculating the specific Lambda components for the PE3
parameter estimation.

## Usage

``` r
cal_Lam(alpha = NULL, eta = 1, r = 1)
```

## Arguments

- alpha:

  A numeric shape parameter.

- eta:

  A non-negative integer representing the order of the LH-moments.
  Default is 1.

- r:

  An integer representing the moment order index. Default is 1.

## Value

A numeric value representing the Lambda component.
