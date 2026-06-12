# Calculate Theoretical LH-moments for the Four-Parameter Kappa Distribution

This function computes the theoretical LH-moments and LH-moment ratios
for the four-parameter Kappa (K4D) distribution given its parameters.
The computations are analytically derived based on the formulas
presented by Murshed et al. (2015).

## Usage

``` r
lhmom.kap(para = NULL, eta = 1)
```

## Arguments

- para:

  A numeric vector of four parameters: c(xi, alpha, k, h), corresponding
  to the location, scale, and the two shape parameters of the K4D
  distribution.

- eta:

  A non-negative integer (between 0 and 4) representing the order of the
  LH-moments. Default is 1.

## Value

A list containing:

- `lambdas`: A named numeric vector of the first four theoretical
  LH-moments.

- `ratios`: A named numeric vector of the corresponding LH-moment
  ratios.

- `eta`: The order of the LH-moments used.

- `type`: The distribution type (`"kap"`).

## References

Murshed, S. M., et al. (2015). Theoretical LH-moments of the Kappa
distribution. *Stochastic Environmental Research and Risk Assessment
(SERRA)*.
