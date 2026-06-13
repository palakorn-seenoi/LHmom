# Estimate Parameters of the Generalized Gumbel (GGD) Distribution using LH-moments

This function estimates the parameters of the Generalized Gumbel (GGD)
distribution based on the sample LH-moments. It evaluates the GGD as a
special case of the three-parameter Kappa distribution where the shape
parameter `k` is fixed at 0.

## Usage

``` r
lh.parggd(data, eta = 1)
```

## Arguments

- data:

  A numeric vector of data values.

- eta:

  A non-negative integer (between 0 and 4) representing the order of the
  LH-moments. Default is 1.

## Value

A list containing:

- `para`: A named numeric vector containing the estimated parameters
  (`mu` for location, `sigma` for scale, and `h` for shape).

- `type`: The distribution type (`"ggd"`).

- `source`: The name of the function (`"lh.parggd"`).

- `eta`: The order of the LH-moments used.
