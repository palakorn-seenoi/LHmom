# Estimate Parameters of the Generalized Pareto (GPA) Distribution using LH-moments

This function estimates the parameters of the Generalized Pareto (GPA)
distribution based on the sample LH-moments. The analytical formulas for
the parameter estimation are derived based on the methodology presented
by Meshgi and Khalili (2009).

## Usage

``` r
lh.pargpa(data, eta = 1)
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
  (`mu` for location, `sigma` for scale, `k` for shape).

- `eta`: The order of the LH-moments used.

- `type`: The distribution type (`"gpa"`).

- `ifail`: A numeric indicator of success (0 for success).

- `source`: The name of the function (`"lh.pargpa"`).

## References

Meshgi, A., & Khalili, D. (2009). Comprehensive evaluation of regional
flood frequency analysis by L- and LH-moments. *Stochastic Environmental
Research and Risk Assessment (SERRA)*, 23(1), 137-152.
