# Estimate Parameters of the Three-Parameter Kappa Distribution with Fixed k

This function estimates the parameters of the three-parameter Kappa
distribution based on the sample LH-moments, given a fixed value for the
shape parameter `k`. It utilizes numerical optimization (`nleqslv`) to
solve for the second shape parameter `h`. If the numerical solver fails
to converge, the function implements a fallback mechanism by adopting
the `h` parameter estimated from the four-parameter Kappa distribution
(`lh.parkap`).

## Usage

``` r
lh.park3d.kfix(data, eta = 1, kfix = 0, hlow = NULL, ntry = 10)
```

## Arguments

- data:

  A numeric vector of data values.

- eta:

  A non-negative integer (between 0 and 4) representing the order of the
  LH-moments. Default is 1.

- kfix:

  A numeric scalar representing the fixed value for the shape parameter
  `k`. Default is 0.

- hlow:

  A numeric scalar representing the lower bound for the shape parameter
  `h`. If `NULL`, it defaults to `-eta - 1`.

- ntry:

  An integer specifying the maximum number of initialization attempts
  for the numerical optimization solver. Default is 10.

## Value

A list containing:

- `para`: A named numeric vector containing the estimated parameters
  (`mu` for location, `sigma` for scale, `kfix`, and `h` for shape).

- `eta`: The order of the LH-moments used.

- `kfix`: The fixed value of the shape parameter `k` used in the
  estimation.

- `type`: The distribution type (`"kap"`).

- `ifail`: A numeric indicator of the optimization solver's status (0
  for success, 2 for partial success/limit reached, 5 for failure).

- `precision`: The final function value (`fvec`) from the solver,
  indicating the precision of the root found.

- `source`: The name of the function (`"lh.park3d.kfix"`).
