# Theoretical LH-moments for the Three-Parameter Kappa Distribution with Fixed h

This function computes the theoretical LH-moments for the
three-parameter Kappa distribution where the shape parameter h is fixed.

## Usage

``` r
lhmom.k3d.hfix(para = NULL, eta = 1, hfix = 0)
```

## Arguments

- para:

  A numeric vector of three parameters: c(mu, sigma, k).

- eta:

  The order of LH-moments (default is 1).

- hfix:

  The fixed numeric value for the shape parameter h (default is 0).

## Value

A list containing the calculated LH-moments and the distribution type
("kap").
