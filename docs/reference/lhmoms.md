# Calculate Sample LH-moments

This function computes the sample LH-moments and LH-moment ratios from a
numeric dataset. The estimation is based on the methodology proposed by
Wang (1997). The order parameter `eta` determines the weight assigned to
larger observations. When `eta = 0`, the calculation reduces to the
ordinary sample L-moments.

## Usage

``` r
lhmoms(data, eta = NULL, nmom = 5)
```

## Arguments

- data:

  A numeric vector of data values.

- eta:

  A non-negative integer (or a sequence of integers) between 0 and 4
  representing the order of the LH-moments.

- nmom:

  An integer specifying the maximum number of moments to compute
  (default is 5).

## Value

A list containing:

- `eta`: The order(s) of the LH-moments.

- `lambdas`: A matrix of the calculated sample LH-moments.

- `ratios`: A matrix of the calculated sample LH-moment ratios.

## References

Wang, Q. J. (1997). Using higher order L-moments for regional flood
frequency analysis. *Water Resources Research*, 33(12), 2841-2848.
