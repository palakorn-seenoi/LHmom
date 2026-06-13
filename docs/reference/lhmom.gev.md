# Calculate Theoretical LH-moments for the Generalized Extreme Value (GEV) Distribution

This function computes the theoretical LH-moments and LH-moment ratios
for the Generalized Extreme Value (GEV) distribution given its
parameters. The computations are derived based on the formulas provided
by Wang (1997).

## Usage

``` r
lhmom.gev(para, eta = NULL)
```

## Arguments

- para:

  A numeric vector of three parameters: c(xi, alpha, k), corresponding
  to the location, scale, and shape parameters of the GEV distribution,
  respectively.

- eta:

  A non-negative integer (between 0 and 4) representing the order of the
  LH-moments.

## Value

A list containing:

- `lambdas`: A named numeric vector of the first four theoretical
  LH-moments.

- `ratios`: A named numeric vector of the corresponding LH-moment
  ratios.

- `eta`: The order of the LH-moments used.

## References

Wang, Q. J. (1997). Using higher order L-moments for regional flood
frequency analysis. *Water Resources Research*, 33(12), 2841-2848.
