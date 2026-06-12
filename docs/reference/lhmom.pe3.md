# Calculate Theoretical LH-moments for the Pearson Type III (PE3) Distribution

This function computes the theoretical LH-moments and LH-moment ratios
for the Pearson Type III (PE3) distribution. It supports parameter
inputs either as a `vec2par` object from the `lmomco` package or as a
standard numeric vector.

## Usage

``` r
lhmom.pe3(para, eta = 0)
```

## Arguments

- para:

  A list containing a `vec2par` object (with `mu, sigma, gamma`) OR a
  numeric vector of length 3 `c(xi, alpha, beta)`.

- eta:

  A non-negative integer representing the order of the LH-moments.
  Default is 0 (which corresponds to the standard L-moments).

## Value

A list containing:

- `lambdas`: A named numeric vector of the first four theoretical
  LH-moments.

- `ratios`: A named numeric vector of the corresponding LH-moment
  ratios.

- `trim`: The trim level (fixed at 0 for LH-moments compatibility).

- `type`: The distribution type (`"pe3"`).

- `eta`: The order of the LH-moments used.

- `source`: The name of the function (`"lhmom.pe3"`).
