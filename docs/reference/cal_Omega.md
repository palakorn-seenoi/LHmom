# Calculate the Omega Function for the PE3 Distribution

An internal helper function to calculate the generalized Probability
Weighted Moments (PWM) component, Omega, for the Pearson Type III (PE3)
distribution via numerical integration.

## Usage

``` r
cal_Omega(alpha = NULL, S = 1)
```

## Arguments

- alpha:

  A numeric shape parameter for the Gamma distribution component.

- S:

  A numeric value representing the moment order step (`eta + k`).
  Default is 1.

## Value

A numeric value resulting from the numerical integration.
