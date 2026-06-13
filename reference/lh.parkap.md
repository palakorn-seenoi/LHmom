# Estimate Parameters of the Four-Parameter Kappa Distribution using LH-moments

This function estimates the parameters of the four-parameter Kappa (K4D)
distribution based on the sample LH-moments. It utilizes numerical
optimization (`nleqslv`) to simultaneously solve for the two shape
parameters, `k` and `h`. If the optimization fails or becomes unstable,
the function provides fallback mechanisms utilizing fixed-parameter
Kappa estimations. When `eta = 0`, it defaults to the ordinary L-moments
estimation via
[`lmomco::parkap`](https://rdrr.io/pkg/lmomco/man/parkap.html).

## Usage

``` r
lh.parkap(
  data,
  eta = 1,
  snap.tau4 = TRUE,
  nudge.tau4 = 1e-05,
  hlow = NULL,
  ntry = 10
)
```

## Arguments

- data:

  A numeric vector of data values.

- eta:

  A non-negative integer (between 0 and 4) representing the order of the
  LH-moments. Default is 1.

- snap.tau4:

  A logical value indicating whether to snap the sample L-kurtosis
  (`tau4`) downward if it slightly exceeds the theoretical upper bound.
  Default is `TRUE`.

- nudge.tau4:

  A small numeric value used to adjust `tau4` downward if
  `snap.tau4 = TRUE`. Default is `1e-5`.

- hlow:

  A numeric scalar representing the lower bound for the shape parameter
  `h`. If `NULL`, it defaults to `-eta - 1`.

- ntry:

  An integer specifying the maximum number of initialization attempts
  for the numerical optimization solver. Default is 10.

## Value

A list containing:

- `type`: The distribution type (`"kap"`).

- `para`: A named numeric vector containing the estimated parameters
  (`mu` for location, `sigma` for scale, `k` and `h` for shapes).

- `eta`: The order of the LH-moments used.

- `source`: The name of the function (`"lh.parkap"`).

- `ifail`: A numeric indicator of the optimization solver's status (0
  for success, 1 for fallback success, 5 for failure).

- `precision`: The final function values from the solver.

- `ifailtext`: A descriptive message regarding the estimation success or
  failure.
