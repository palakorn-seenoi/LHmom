# Theoretical LH-moments for Gumbel Distribution

This function computes the theoretical LH-moments for the Gumbel
distribution by evaluating it as a special case of the Generalized
Extreme Value (GEV) distribution with the shape parameter set to zero.

## Usage

``` r
lhmom.gum(para = NULL, eta = 1)
```

## Arguments

- para:

  A vector of parameters c(mu, sigma)

- eta:

  The order of LH-moments (default is 1).

## Value

A list containing the LH-moments and distribution type.
