# Estimate Parameters of the Generalized Extreme Value (GEV) Distribution using LH-moments

This function estimates the parameters of the Generalized Extreme Value
(GEV) distribution based on the sample LH-moments. The estimation
methodology follows Wang (1997). It provides two approaches for
estimating the shape parameter: using Wang's predefined polynomial
approximations or utilizing numerical optimization.

## Usage

``` r
lh.pargev(data, eta = 1, opt = FALSE, ntry = 5)
```

## Arguments

- data:

  A numeric vector of data values.

- eta:

  A non-negative integer (between 0 and 4) representing the order of the
  LH-moments. Default is 1.

- opt:

  A logical value indicating the estimation method for the shape
  parameter. If `FALSE` (default), it uses Wang's polynomial
  approximation coefficients. If `TRUE`, it uses numerical optimization
  via `nleqslv`.

- ntry:

  An integer specifying the maximum number of initialization attempts
  for the numerical optimization solver when `opt = TRUE`. Default is 5.

## Value

A list containing:

- `type`: The distribution type (`"gev"`).

- `para`: A named numeric vector containing the estimated parameters
  (`xi` for location, `alpha` for scale, `k` for shape).

- `eta`: The order of the LH-moments used.

- `source`: The name of the function (`"lh.pargev"`).

- `ifail`: A numeric indicator of the optimization solver's status (0
  for success, 5 for failure to converge).

## References

Wang, Q. J. (1997). Using higher order L-moments for regional flood
frequency analysis. *Water Resources Research*, 33(12), 2841-2848.
