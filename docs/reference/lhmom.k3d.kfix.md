# Theoretical LH-moments for the Three-Parameter Kappa Distribution with Fixed k

This function computes the theoretical LH-moments for the
three-parameter Kappa distribution where the shape parameter k is fixed.

## Usage

``` r
lhmom.k3d.kfix(para = NULL, eta = 1, kfix = 0)
```

## Arguments

- para:

  A numeric vector of three parameters: c(mu, sigma, h).

- eta:

  The order of LH-moments (default is 1).

- kfix:

  The fixed numeric value for the shape parameter k (default is 0).

## Value

A list containing the calculated LH-moments and the distribution type
("kap").
