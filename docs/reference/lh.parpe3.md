# Estimate Parameters of the Pearson Type III (PE3) Distribution using LH-moments

This function estimates the parameters of the Pearson Type III (PE3)
distribution based on the sample LH-moments. It provides two estimation
methods: using predefined polynomial coefficients (Log-Log Split Degree
5 based on sample LH-skewness) or numerical optimization via `nleqslv`.
If the numerical solver fails to converge, the function falls back to
the ordinary L-moments estimation (`eta = 0`).

## Usage

``` r
lh.parpe3(data, eta = 1, opt = FALSE, ntry = 10)
```

## Arguments

- data:

  A numeric vector of data values.

- eta:

  A non-negative integer (between 0 and 4) representing the order of the
  LH-moments. Default is 1.

- opt:

  A logical value. If `FALSE` (default), it estimates parameters using
  polynomial approximations. If `TRUE`, it utilizes numerical
  optimization.

- ntry:

  An integer specifying the maximum number of initialization attempts
  for the numerical optimization solver when `opt = TRUE`. Default is
  10.

## Value

A list containing:

- `type`: The distribution type (`"pe3"`).

- `para`: A named numeric vector containing the estimated parameters
  (`mu` for location, `sigma` for scale, `gamma` for shape).

- `eta`: The order of the LH-moments used.

- `source`: The name of the function (`"lh.parpe3"`).

- `ifail`: A numeric indicator of the optimization solver's status (0
  for success, 2 for partial success, 5 for failure).

- `precision`: The final function value (`fvec`) from the solver.
