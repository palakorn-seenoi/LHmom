# Estimate Parameters of the Three-Parameter Kappa Distribution with Fixed h

This function estimates the parameters of the three-parameter Kappa
distribution based on the sample LH-moments, given a fixed value for the
shape parameter `h`. The function utilizes numerical optimization
(`nleqslv`) to solve for the shape parameter `k`. If the numerical
solver fails to converge, the function implements a fallback mechanism
by adopting the `k` parameter estimated from the Generalized Extreme
Value (GEV) distribution.

## Usage

``` r
lh.park3d.hfix(data, eta = 1, hfix = 0, ntry = 10)
```

## Arguments

- data:

  A numeric vector of data values.

- eta:

  A non-negative integer (between 0 and 4) representing the order of the
  LH-moments. Default is 1.

- hfix:

  A numeric scalar representing the fixed value for the shape parameter
  `h`. Default is 0.

- ntry:

  An integer specifying the maximum number of initialization attempts
  for the numerical optimization solver. Default is 10.

## Value

A list containing:

- `para`: A named numeric vector containing the estimated parameters
  (`mu` for location, `sigma` for scale, `k` for shape, and `hfix`).

- `eta`: The order of the LH-moments used.

- `hfix`: The fixed value of the shape parameter `h` used in the
  estimation.

- `type`: The distribution type (`"kap"`).

- `ifail`: A numeric indicator of the optimization solver's status (0
  for success, 2 for partial success/limit reached, 5 for failure).

- `precision`: The final function value (`fvec`) from the solver,
  indicating the precision of the root found.

- `source`: The name of the function (`"lh.park3d.hfix"`).
