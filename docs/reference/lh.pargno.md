# Estimate Parameters of the Generalized Normal (GNO) Distribution using LH-moments

This function estimates the parameters of the Generalized Normal (GNO)
distribution based on the sample LH-moments. It provides two methods for
estimating the shape parameter: using predefined polynomial
approximations or utilizing numerical optimization.

## Usage

``` r
lh.pargno(data, eta = 1, opt = FALSE, ntry = 5)
```

## Arguments

- data:

  A numeric vector of data values.

- eta:

  A non-negative integer (between 0 and 4) representing the order of the
  LH-moments. Default is 1.

- opt:

  A logical value indicating the estimation method for the shape
  parameter. If `FALSE` (default), it uses a polynomial approximation.
  If `TRUE`, it uses numerical optimization via `nleqslv`.

- ntry:

  An integer specifying the maximum number of initialization attempts
  for the numerical optimization solver when `opt = TRUE`. Default is 5.

## Value

A list containing:

- `type`: The distribution type (`"gno"`).

- `para`: A named numeric vector containing the estimated parameters
  (`xi` for location, `alpha` for scale, `k` for shape).

- `eta`: The order of the LH-moments used.

- `source`: The name of the function (`"lh.pargno"`).

- `ifail`: A numeric indicator of the optimization solver's status (0
  for success, 5 for failure to converge).
