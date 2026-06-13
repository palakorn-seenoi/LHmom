# Generate Initial Parameters for Four-Parameter Kappa Optimization

An internal helper function used to initialize the shape parameters (`k`
and `h`) for the numerical optimization routines involved in estimating
the four-parameter Kappa distribution via LH-moments.

## Usage

``` r
initkh.k4d(data, ntry = 5)
```

## Arguments

- data:

  A numeric vector of data values.

- ntry:

  An integer specifying the number of initial parameter combinations to
  generate. Minimum is 5. Default is 5.

## Value

A matrix containing `ntry` rows and 2 columns, representing various
starting points for the `k` and `h` shape parameters.
