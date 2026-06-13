# Theoretical LH-moments for the Generalized Pareto (GPA) Distribution

This function computes the theoretical LH-moments for the Generalized
Pareto (GPA) distribution. It is evaluated as a special case of the
four-parameter Kappa distribution with the shape parameter h set to 1.

## Usage

``` r
lhmom.gpa(para = NULL, eta = 1)
```

## Arguments

- para:

  A numeric vector of three parameters: c(mu, sigma, k).

- eta:

  The order of LH-moments (default is 1).

## Value

A list containing the calculated LH-moments and the distribution type
("gpa").
