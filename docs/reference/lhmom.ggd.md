# Theoretical LH-moments for the Generalized Gumbel (GGD) Distribution

This function computes the theoretical LH-moments for the Generalized
Gumbel (GGD) distribution. It is evaluated as a special case of the
four-parameter Kappa distribution with the shape parameter k set to 0.

## Usage

``` r
lhmom.ggd(para = NULL, eta = 1)
```

## Arguments

- para:

  A numeric vector of three parameters: c(mu, sigma, h).

- eta:

  The order of LH-moments (default is 1).

## Value

A list containing the calculated LH-moments and the distribution type
("ggd").
