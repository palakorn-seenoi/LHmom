# Estimate Parameters of the Gumbel Distribution using LH-moments

This function estimates the parameters of the Gumbel distribution based
on the sample LH-moments. The Gumbel distribution is treated as a
special case of the Generalized Extreme Value (GEV) distribution with
the shape parameter (k) fixed at 0. The analytical estimation utilizes
the first two LH-moments and the Euler-Mascheroni constant.

## Usage

``` r
lh.pargum(data, eta = 1)
```

## Arguments

- data:

  A numeric vector of data values.

- eta:

  A non-negative integer (between 0 and 4) representing the order of the
  LH-moments. Default is 1.

## Value

A list containing:

- `type`: The distribution type (`"gum"`).

- `para`: A named numeric vector containing the estimated parameters
  (`xi` for location, `alpha` for scale, and `k` = 0 for shape).

- `eta`: The order of the LH-moments used.

- `source`: The name of the function (`"lh.pargum"`).

- `ifail`: A numeric indicator of success (0 for success).
