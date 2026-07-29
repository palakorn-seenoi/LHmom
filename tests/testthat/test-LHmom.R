library(testthat)
library(LHmom)



# =====================================================================
# 1. Test: Lock in the sample LH-moment values for bundled datasets
# =====================================================================

test_that("Sample LH-moments are calculated correctly for bangkok1 (eta = 1)", {
  data(bangkok1)
  sample_lh <- lhmoms(bangkok1$rainfall, eta = 1)

  expect_equal(as.numeric(sample_lh$lambdas[1]), 126.5447, tolerance = 1e-4)
  expect_equal(as.numeric(sample_lh$lambdas[2]), 18.86068, tolerance = 1e-4)
})

test_that("Sample LH-moments are calculated correctly for khonkaen (eta = 1)", {
  data(khonkaen)
  sample_lh <- lhmoms(khonkaen$rainfall, eta = 1)

  expect_equal(as.numeric(sample_lh$lambdas[1]), 117.35923, tolerance = 1e-4)
  expect_equal(as.numeric(sample_lh$lambdas[2]), 19.5594,   tolerance = 1e-4)
})

test_that("Sample LH-moments are calculated correctly for sarakham (eta = 1)", {
  data(sarakham)
  sample_lh <- lhmoms(sarakham$temperature, eta = 1)

  expect_equal(as.numeric(sample_lh$lambdas[1]), 41.3937805, tolerance = 1e-4)
  expect_equal(as.numeric(sample_lh$lambdas[2]), 0.4505488,  tolerance = 1e-4)
})


test_that("Sample LH-moments are calculated correctly for bangkok1 (eta = 2)", {
  data(bangkok1)
  sample_lh <- lhmoms(bangkok1$rainfall, eta = 2)

  expect_equal(as.numeric(sample_lh$lambdas[1]), 139.11847, tolerance = 1e-4)
  expect_equal(as.numeric(sample_lh$lambdas[2]), 18.22076,  tolerance = 1e-4)
})

test_that("Sample LH-moments are calculated correctly for khonkaen (eta = 2)", {
  data(khonkaen)
  sample_lh <- lhmoms(khonkaen$rainfall, eta = 2)

  expect_equal(as.numeric(sample_lh$lambdas[1]), 130.39889, tolerance = 1e-4)
  expect_equal(as.numeric(sample_lh$lambdas[2]), 19.84207,  tolerance = 1e-4)
})

test_that("Sample LH-moments are calculated correctly for sarakham (eta = 2)", {
  data(sarakham)
  sample_lh <- lhmoms(sarakham$temperature, eta = 2)

  expect_equal(as.numeric(sample_lh$lambdas[1]), 41.6941463, tolerance = 1e-4)
  expect_equal(as.numeric(sample_lh$lambdas[2]), 0.4018386,  tolerance = 1e-4)
})

# =====================================================================
# 2. Test: Round-trip consistency for 3-Parameter Models (GEV)
# =====================================================================

test_that("Round-trip consistency holds for GEV distribution for bangkok1 at eta = 1", {
  data(bangkok1)
  sample_lh <- lhmoms(bangkok1$rainfall,    eta = 1)
  fit_gev   <- lh.pargev(bangkok1$rainfall, eta = 1)
  theo_lh   <- lhmom.gev(fit_gev$para,      eta = 1)

  # First two raw moments must match perfectly
  expect_equal(as.numeric(theo_lh$lambdas[1:2]), as.numeric(sample_lh$lambdas[1:2]), tolerance = 1e-6)

  # Third ratio must be close (approximation)
  expect_equal(as.numeric(theo_lh$ratios[3]), as.numeric(sample_lh$ratios[3]), tolerance = 5e-3)
})

test_that("Round-trip consistency holds for GEV distribution for khonkaen at eta = 1", {
  data(khonkaen)
  sample_lh <- lhmoms(khonkaen$rainfall,    eta = 1)
  fit_gev   <- lh.pargev(khonkaen$rainfall, eta = 1)
  theo_lh   <- lhmom.gev(fit_gev$para,      eta = 1)

  # First two raw moments must match perfectly
  expect_equal(as.numeric(theo_lh$lambdas[1:2]), as.numeric(sample_lh$lambdas[1:2]), tolerance = 1e-6)

  # Third ratio must be close (approximation)
  expect_equal(as.numeric(theo_lh$ratios[3]), as.numeric(sample_lh$ratios[3]), tolerance = 5e-3)
})

test_that("Round-trip consistency holds for GEV distribution for sarakham at eta = 1", {
  data("sarakham")
  sample_lh <- lhmoms(sarakham$temperature,    eta = 1)
  fit_gev   <- lh.pargev(sarakham$temperature, eta = 1)
  theo_lh   <- lhmom.gev(fit_gev$para,         eta = 1)

  # First two raw moments must match perfectly
  expect_equal(as.numeric(theo_lh$lambdas[1:2]), as.numeric(sample_lh$lambdas[1:2]), tolerance = 1e-6)

  # Third ratio must be close (approximation)
  expect_equal(as.numeric(theo_lh$ratios[3]), as.numeric(sample_lh$ratios[3]), tolerance = 5e-3)
})




# =====================================================================
# 3. Test: Round-trip consistency for 4-Parameter Models (Kappa)
# =====================================================================

test_that("Round-trip consistency holds for Kappa distribution for bangkok1 at eta = 1", {
  data(bangkok1)
  sample_lh <- lhmoms(bangkok1$rainfall,    eta = 1)
  fit_kap <-   lh.parkap(bangkok1$rainfall, eta = 1)
  theo_lh <-   lhmom.kap(fit_kap$para,      eta = 1)

  # First two raw moments and subsequent two ratios must match
  expect_equal(as.numeric(theo_lh$lambdas[1:2]), as.numeric(sample_lh$lambdas[1:2]), tolerance = 1e-6)
  expect_equal(as.numeric(theo_lh$ratios[3:4]),  as.numeric(sample_lh$ratios[3:4]),  tolerance = 1e-4)
})


test_that("Round-trip consistency holds for Kappa distribution for khonkaen at eta = 1", {
  data(khonkaen)
  sample_lh <- lhmoms(khonkaen$rainfall,    eta = 1)
  fit_kap <-   lh.parkap(khonkaen$rainfall, eta = 1)
  theo_lh <-   lhmom.kap(fit_kap$para,      eta = 1)

  # First two raw moments and subsequent two ratios must match
  expect_equal(as.numeric(theo_lh$lambdas[1:2]), as.numeric(sample_lh$lambdas[1:2]), tolerance = 1e-6)
  expect_equal(as.numeric(theo_lh$ratios[3:4]),  as.numeric(sample_lh$ratios[3:4]),  tolerance = 1e-4)
})


test_that("Round-trip consistency holds for Kappa distribution for sarakham at eta = 1", {
  data(sarakham)
  sample_lh <- lhmoms(sarakham$temperature,    eta = 1)
  fit_kap <-   lh.parkap(sarakham$temperature, eta = 1)
  theo_lh <-   lhmom.kap(fit_kap$para,         eta = 1)

  # First two raw moments and subsequent two ratios must match
  expect_equal(as.numeric(theo_lh$lambdas[1:2]), as.numeric(sample_lh$lambdas[1:2]), tolerance = 1e-6)
  expect_equal(as.numeric(theo_lh$ratios[3:4]),  as.numeric(sample_lh$ratios[3:4]),  tolerance = 1e-4)
})


test_that("Round-trip consistency holds for Kappa distribution for bangkok1 at eta = 2", {
  data(bangkok1)
  sample_lh <- lhmoms(bangkok1$rainfall,    eta = 2)
  fit_kap <-   lh.parkap(bangkok1$rainfall, eta = 2)
  theo_lh <-   lhmom.kap(fit_kap$para,      eta = 2)

  # First two raw moments and subsequent two ratios must match
  expect_equal(as.numeric(theo_lh$lambdas[1:2]), as.numeric(sample_lh$lambdas[1:2]), tolerance = 1e-6)
  expect_equal(as.numeric(theo_lh$ratios[3:4]),  as.numeric(sample_lh$ratios[3:4]),  tolerance = 1e-4)
})


test_that("Round-trip consistency holds for Kappa distribution for khonkaen at eta = 2", {
  data(khonkaen)
  sample_lh <- lhmoms(khonkaen$rainfall,    eta = 2)
  fit_kap <-   lh.parkap(khonkaen$rainfall, eta = 2)
  theo_lh <-   lhmom.kap(fit_kap$para,      eta = 2)

  # First two raw moments and subsequent two ratios must match
  expect_equal(as.numeric(theo_lh$lambdas[1:2]), as.numeric(sample_lh$lambdas[1:2]), tolerance = 1e-6)
  expect_equal(as.numeric(theo_lh$ratios[3:4]),  as.numeric(sample_lh$ratios[3:4]),  tolerance = 1e-4)
})


test_that("Round-trip consistency holds for Kappa distribution for sarakham at eta = 2", {
  data(sarakham)
  sample_lh <- lhmoms(sarakham$temperature,    eta = 2)
  fit_kap <-   lh.parkap(sarakham$temperature, eta = 2)
  theo_lh <-   lhmom.kap(fit_kap$para,         eta = 2)

  # First two raw moments and subsequent two ratios must match
  expect_equal(as.numeric(theo_lh$lambdas[1:2]), as.numeric(sample_lh$lambdas[1:2]), tolerance = 1e-6)
  expect_equal(as.numeric(theo_lh$ratios[3:4]),  as.numeric(sample_lh$ratios[3:4]),  tolerance = 1e-4)
})


