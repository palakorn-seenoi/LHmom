# Quantile-Quantile Plot for LH-moments Fitted Distributions

Generates a Quantile-Quantile (Q-Q) plot to visually assess the
goodness-of-fit of a probability distribution fitted using higher-order
L-moments (LH-moments). The function plots the empirical quantiles of
the data against the theoretical quantiles of the specified
distribution. It optionally calculates and displays bootstrap confidence
intervals to evaluate model uncertainty.

## Usage

``` r
lh.qqplot(
  data,
  fit_obj,
  main = "Q-Q Plot",
  ci = FALSE,
  ci.level = 0.95,
  n.boot = 500
)
```

## Arguments

- data:

  A numeric vector of observations or a data frame with 1 column.

- fit_obj:

  An object returned by any `lh.par*` function containing the fitted
  parameters.

- main:

  Title of the plot.

- ci:

  Logical; if `TRUE`, calculates and plots the bootstrap confidence
  interval.

- ci.level:

  Numeric; confidence level for the interval (e.g., 0.95 for 95%).

- n.boot:

  Number of bootstrap samples for CI (default is 500).

## Value

An invisible data frame containing theoretical quantiles, empirical
quantiles, and CI bounds (if calculated).
