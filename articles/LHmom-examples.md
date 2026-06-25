# Illustrative Examples of LHmom

## Introduction

This vignette demonstrates the core functionalities of the `LHmom`
package using the built-in `bangkok` dataset. This dataset contains the
annual maximum series of daily rainfall in Bangkok, Thailand, from 1951
to 2025.

We will walk through the process of calculating sample LH-moments,
estimating parameters for the Generalized Extreme Value (GEV)
distribution, and evaluating the model fit.

## 1. Installation

You can install the development version of LHmom from GitHub with:

``` r

install.packages("devtools")
devtools::install_github("palakorn-seenoi/LHmom")
```

## 2. Data Preparation

Load the package and the dataset:

``` r

library(LHmom)
data(bangkok1)

# View the first few rows
head(bangkok1)
#>   year rainfall
#> 1 1951    133.5
#> 2 1952    111.0
#> 3 1953     84.1
#> 4 1954     53.8
#> 5 1955    108.8
#> 6 1956     69.4
```

## 3. LH estimation for GEV distribution

Users can compute the sample LH-moments for a given order $`\eta`$ and
subsequently estimate the parameters for various distributions. The
following example demonstrates fitting the GEV distribution using the
first order of LH-moments ($`\eta = 1`$).

### Calculating Sample LH-moments

``` r

# Calculate sample LH-moments at eta = 1
sample_lh <- lhmoms(bangkok1$rainfall, eta = 1)
print(sample_lh$lambdas)
#>        lhmom-1  lhmom-2  lhmom-3  lhmom-4  lhmom-5
#> eta=1 126.5447 18.86068 5.220361 3.037946 1.666922
```

### Parameter Estimation

Next, we estimate the GEV parameters using the calculated sample
LH-moments:

``` r

# Estimate GEV parameters (eta = 1)
fit_gev <- lh.pargev(bangkok1$rainfall, eta = 1)
print(fit_gev$para)
#>          xi       alpha           k 
#> 87.64384852 28.23812277 -0.06139394
```

### Verifying Numerical Consistency

To verify the estimation algorithm, we can calculate the theoretical
LH-moments from our newly fitted parameters:

``` r

# Calculate theoretical LH-moments from the fitted parameters
theo_lh <- lhmom.gev(fit_gev$para, eta = 1)
print(theo_lh$lambdas)
#>    LHmom-1    LHmom-2    LHmom-3    LHmom-4 
#> 126.544685  18.860680   5.229035   3.131879
```

Notice the numerical consistency: because the GEV is a three-parameter
distribution, the estimation algorithm successfully equates the first
three sample LH-moments with their theoretical counterparts.

## 4. Wang’s Goodness-of-Fit Test

A crucial aspect of extreme value analysis is model validation. The
[`wang.test.lhgev()`](https://palakorn-seenoi.github.io/LHmom/reference/wang.test.lhgev.md)
function evaluates whether the GEV distribution provides a suitable fit
to the data across different LH-moment orders ($`\eta`$) by calculating
a Z-statistic based on the difference between the sample and theoretical
LH-kurtosis.

``` r

# Perform Wang's goodness-of-fit test
gof_test <- wang.test.lhgev(bangkok1$rainfall)
print(gof_test)
#>   eta       z.test cond.sigma   p.value
#> 1   0 -0.149710710 0.03581892 0.8809929
#> 2   1 -0.167327234 0.02976428 0.8671126
#> 3   2 -0.170480325 0.02714540 0.8646324
#> 4   3  0.006399968 0.02572831 0.9948936
#> 5   4  0.170936924 0.02478531 0.8642734
```

### Interpretation:

The results demonstrate a satisfactory fit across all $`\eta`$ values,
as all $`p`$-values significantly exceed the 0.05 significance level.
Notably, at $`\eta = 3`$, the absolute Z-statistic reaches its minimum,
and the $`p`$-value approaches unity. This indicates that employing a
higher-order LH-moment ($`\eta = 3`$) captures the upper-tail
distributional shape of the annual maximum rainfall exceptionally well
compared to ordinary L-moments ($`\eta = 0`$).
