# Wang's Goodness-of-Fit Test for the Generalized Extreme Value (GEV) Distribution

This function performs a Goodness-of-Fit (GOF) test for the Generalized
Extreme Value (GEV) distribution using LH-moments, based on the
methodology proposed by Wang (1998). It calculates a Z-test statistic by
comparing the sample LH-kurtosis with the theoretical LH-kurtosis. The
test is evaluated across LH-moment orders (`eta`) from 0 to 4.

## Usage

``` r
wang.test.lhgev(data)
```

## Arguments

- data:

  A numeric vector of data values.

## Value

A data frame containing the GOF test results for each `eta` value (from
0 to 4), with the following columns:

- `eta`: The order of the LH-moments.

- `z.test`: The calculated Z-statistic for the GOF test.

- `cond.sigma`: The conditional standard deviation of the sample
  LH-kurtosis.

- `p.value`: The two-sided p-value corresponding to the Z-test
  statistic.

## References

Wang, Q. J. (1998). Approximate goodness-of-fit tests of fitted
generalized extreme value distributions using LH moments. *Water
Resources Research*, 34(12), 3497-3502.
