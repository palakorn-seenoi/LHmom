# Calculate Theoretical LH-moments for the Generalized Normal (GNO) Distribution

This function computes the theoretical LH-moments and LH-moment ratios
for the Generalized Normal (GNO) distribution given its parameters. When
the order `eta = 0`, it falls back to the ordinary L-moments using the
`lmomco` package.

## Usage

``` r
lhmom.gno(para = NULL, eta = 1)
```

## Arguments

- para:

  A numeric vector of three parameters: c(xi, alpha, k), corresponding
  to the location, scale, and shape parameters of the GNO distribution,
  respectively.

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

- `type`: The distribution type (`"gno"`).
