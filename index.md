# LHmom: Parameter Estimation and Calculations of Linear Higher-Order Moments

The **LHmom** package computes sample and theoretical Linear
Higher-Order Moments (LH-moments) and implements parameter estimation
for various statistical distributions, including the Generalized Extreme
Value (GEV), Generalized Logistic (GLO), Generalized Pareto (GPA),
Generalized Normal (GNO), Pearson Type III (PE3), and Kappa
distributions. It also provides analytical tools such as Wang’s
goodness-of-fit test.

The design, parameterization, and output structures of this package
closely follow the style and conventions of the ‘lmomco’ package to
ensure a familiar and consistent interface for users accustomed to
L-moment calculations in R.

------------------------------------------------------------------------

## 📚 Documentation

Explore the full capabilities of the **LHmom** package through our
detailed guides:

- [**Function
  Reference**](https://palakorn-seenoi.github.io/LHmom/reference/index.md):
  Comprehensive details on all functions, grouped by parameter
  estimation, theoretical calculations, and data sets.
- [**Illustrative Examples
  (Articles)**](https://palakorn-seenoi.github.io/LHmom/articles/index.md):
  A step-by-step tutorial demonstrating how to calculate sample
  LH-moments, estimate GEV parameters, and perform goodness-of-fit tests
  using real data.

------------------------------------------------------------------------

## 📚 Developers & Affiliations

- **Palakorn Seenoi**  
  Department of Statistics, Faculty of Science, Khon Kaen University,
  Thailand  
  Email: <palakorns@kku.ac.th>

- **Yire Shin**  
  Department of Statistics, Chonnam National University, South Korea  
  Email: <shinyire@daum.net>

- **Piyapatr Busababodhin**  
  Department of Mathematics, Faculty of Science, Mahasarakham
  University, Thailand  
  Email: <piyapatr.b@msu.ac.th>

- **Jeong-Soo Park**  
  Department of Statistics, Chonnam National University, South Korea  
  Email: <jspark@jnu.ac.kr>

------------------------------------------------------------------------

## 📚 Installation

You can install the development version of LHmom from GitHub with:

    # install.packages("devtools")
    devtools::install_github("palakorn-seenoi/LHmom")

## 📚 Quick Example

Here is a quick example of how to calculate sample LH-moments and
estimate GEV parameters using the built-in Bangkok rainfall dataset:

    library(LHmom)
    data(bangkok)

    # 1. Compute sample LH-moments (eta = 1)
    sample_lh <- lhmoms(bangkok$rainfall, eta = 1)
    print(sample_lh$lambdas)

    # 2. Estimate GEV parameters
    fit_gev <- lh.pargev(bangkok$rainfall, eta = 1)
    print(fit_gev$para)
