# LHmom: An R Package for Higher-Order L-Moment Estimation

<!-- badges: start -->
[![GitHub release (latest by date)](https://img.shields.io/github/v/release/palakorn-seenoi/LHmom)](https://github.com/palakorn-seenoi/LHmom/releases)
[![License: GPL (>= 3)](https://img.shields.io/badge/License-GPL%20%28%3E%3E%203%29-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

The **LHmom** package provides a comprehensive framework for parameter estimation using higher-order L-moments (LH-moments). 
By assigning greater weight to upper-tail observations, LH-moments offer a robust statistical method for extreme value analysis, particularly in hydrological and atmospheric research.

Furthermore, the design, parameterization, and output structures of this package closely follow the style and conventions of the widely used `lmomco` package. This ensures a familiar, seamless, and consistent interface for users already accustomed to standard L-moment calculations in R.


## Features

* **Parameter Estimation:** Supports LH-moment orders zero through four (`eta = 0` to `4`) for 10 continuous probability distributions, including:
  * Generalized Extreme Value (GEV)
  * Generalized Pareto (GPA)
  * Generalized Logistic (GLO)
  * Generalized Normal (GNO)
  * Pearson Type III (PE3)
  * Gumbel (GUM)
  * Generalized Gumbel (GGD)
  * Three-parameter Kappa (fixed $k$ or $h$)
  * Highly flexible Four-parameter Kappa (KAP)
* **Diagnostic Tools:** Includes universal Quantile-Quantile (Q-Q) plotting for visual goodness-of-fit validation and an explicit implementation of Wang's goodness-of-fit test for the GEV distribution.
* **Built-in Datasets:** Includes historical environmental datasets (`bangkok1`, `khonkaen`, and `sarakham`) for immediate demonstration in flood frequency and extreme climatic event modeling.

## Installation

You can install the development version of `LHmom` from GitHub using the `devtools` package:

```R
# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Install LHmom from GitHub
devtools::install_github("palakorn-seenoi/LHmom")
```

## Quick Start Example
Here is a basic workflow using the built-in `bangkok1` dataset to estimate parameters for the Generalized Extreme Value (GEV) distribution using the first order of LH-moments (`eta = 1`):

```R
library(LHmom)

# Load the built-in annual maximum rainfall dataset
data(bangkok1)

# 1. Calculate sample LH-moments at eta = 1
sample_lh <- lhmoms(bangkok1$rainfall, eta = 1)
print(sample_lh$lambdas)

# 2. Estimate GEV parameters using sample LH-moments
fit_gev <- lh.pargev(bangkok1$rainfall, eta = 1)
print(fit_gev$para)

# 3. Generate Q-Q plot with 95% bootstrap confidence intervals
lh.qqplot(bangkok1$rainfall, fit_gev, 
          main = "GEV Q-Q Plot (eta = 1)", 
          ci = TRUE, ci.level = 0.95)

# 4. Perform Wang's goodness-of-fit test for the GEV distribution
wang.test.lhgev(bangkok1$rainfall)
```

## Documentation

You can find the complete HTML documentation and tutorials at [LHmom Online Manual](https://palakorn-seenoi.github.io/LHmom/).

## Citation

If you use `LHmom` in your research, please cite the following paper:

Seenoi, P., Shin, Y., Busababodhin, P., & Park, J.S. (2026). LHmom: An R Package for higher order L-Moment Estimation. SoftwareX. (Submitted/Under Review)


## Contact
For questions, bug reports, or feature requests, please open an issue on GitHub or contact Palakorn Seenoi at palakorns@kku.ac.th.
