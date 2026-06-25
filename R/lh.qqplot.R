#' Quantile-Quantile Plot for LH-moments Fitted Distributions
#'
#' @description
#' Generates a Quantile-Quantile (Q-Q) plot to visually assess the goodness-of-fit
#' of a probability distribution fitted using higher-order L-moments (LH-moments).
#' The function plots the empirical quantiles of the data against the theoretical
#' quantiles of the specified distribution. It optionally calculates and displays
#' bootstrap confidence intervals to evaluate model uncertainty.
#'
#' @param data A numeric vector of observations or a data frame with 1 column.
#' @param fit_obj An object returned by any \code{lh.par*} function containing the fitted parameters.
#' @param main Title of the plot.
#' @param ci Logical; if \code{TRUE}, calculates and plots the bootstrap confidence interval.
#' @param ci.level Numeric; confidence level for the interval (e.g., 0.95 for 95%).
#' @param n.boot Number of bootstrap samples for CI (default is 500).
#' @param ... Additional graphical parameters passed to \code{plot()}.
#'
#' @importFrom graphics plot abline polygon points legend
#' @importFrom stats na.omit quantile
#' @importFrom grDevices rgb
#'
#' @return An invisible data frame containing theoretical quantiles, empirical quantiles,
#' and CI bounds (if calculated).
#' @export
#'
lh.qqplot <- function(data, fit_obj, main = "Q-Q Plot", ci = FALSE, ci.level = 0.95, n.boot = 500, ...) {

  if (is.data.frame(data) || is.list(data)) data <- unlist(data)
  data <- as.numeric(data)
  data <- sort(na.omit(data))
  n <- length(data)
  p <- (1:n - 0.35) / (n + 0.3)

  dist_type <- tolower(fit_obj$type)
  eta_val <- fit_obj$eta

  get_theo_q <- function(prob, fit) {
    switch(tolower(fit$type),
           "gev"      = lmomco::quagev(prob, fit),
           "glo"      = lmomco::quaglo(prob, fit),
           "gpa"      = lmomco::quagpa(prob, fit),
           "gno"      = lmomco::quagno(prob, fit),
           "pe3"      = lmomco::quape3(prob, fit),
           "gum"      = lmomco::quagum(prob, fit),
           "kap"      = lmomco::quakap(prob, fit),
           "k3d_hfix" = quak3d.hfix(prob, fit),
           "k3d_kfix" = quak3d.kfix(prob, fit),
           "ggd"      = quaggd(prob, fit),
           stop(paste("Error: Distribution type '", fit$type, "' is not supported."))
    )
  }

  theo_q <- get_theo_q(p, fit_obj)

  lower_ci <- rep(NA, n)
  upper_ci <- rep(NA, n)

  if (ci) {
    if (ci.level > 1) ci.level <- ci.level / 100
    if (ci.level <= 0 || ci.level >= 1) stop("Error: ci.level must be between 0 and 1.")

    alpha <- 1 - ci.level
    prob_lower <- alpha / 2
    prob_upper <- 1 - (alpha / 2)

    boot_mat <- matrix(NA, nrow = n, ncol = n.boot)

    for (i in 1:n.boot) {
      tryCatch({
        data_boot <- sample(data, size = n, replace = TRUE)
        fit_boot <- switch(dist_type,
                           "gev"      = lh.pargev(data_boot, eta = eta_val),
                           "glo"      = lh.parglo(data_boot, eta = eta_val),
                           "gpa"      = lh.pargpa(data_boot, eta = eta_val),
                           "gno"      = lh.pargno(data_boot, eta = eta_val),
                           "pe3"      = lh.parpe3(data_boot, eta = eta_val),
                           "gum"      = lh.pargum(data_boot, eta = eta_val),
                           "kap"      = lh.parkap(data_boot, eta = eta_val),
                           "ggd"      = lh.parggd(data_boot, eta = eta_val),
                           "k3d_hfix" = lh.park3d.hfix(data_boot, eta = eta_val, hfix = fit_obj$para[4]),
                           "k3d_kfix" = lh.park3d.kfix(data_boot, eta = eta_val, kfix = fit_obj$para[3]),
                           stop("Refitting not supported.")
        )
        boot_mat[, i] <- get_theo_q(p, fit_boot)
      }, error = function(e) { })
    }
    lower_ci <- apply(boot_mat, 1, quantile, probs = prob_lower, na.rm = TRUE)
    upper_ci <- apply(boot_mat, 1, quantile, probs = prob_upper, na.rm = TRUE)
  }

  x_label <- paste0("Theoretical Quantiles (", toupper(dist_type),
                    ifelse(!is.null(eta_val), paste0(", eta = ", eta_val), ""), ")")

  plot(theo_q, data, type = "n", xlab = x_label, ylab = "Empirical Quantiles", main = main, ...)

  if (ci && !all(is.na(lower_ci))) {
    polygon(c(theo_q, rev(theo_q)), c(lower_ci, rev(upper_ci)), col = rgb(0, 0, 1, 0.15), border = NA)
  }

  abline(0, 1, col = "red", lwd = 2, lty = 2)
  points(theo_q, data, pch = 19, col = "blue")

  # --- Main Legend (Top Left) ---
  if (ci && !all(is.na(lower_ci))) {
    legend("topleft", legend = c("Empirical", "y = x", paste0(round(ci.level * 100), "% CI")),
           col = c("blue", "red", rgb(0, 0, 1, 0.3)), pch = c(19, NA, 15), lty = c(NA, 2, NA),
           lwd = c(NA, 2, NA), pt.cex = c(1, 1, 2), bty = "n", cex = 0.9)
  } else {
    legend("topleft", legend = c("Empirical", "y = x"), col = c("blue", "red"),
           pch = c(19, NA), lty = c(NA, 2), lwd = c(NA, 2), bty = "n", cex = 0.9)
  }

  invisible(data.frame(
    Theoretical = theo_q, Empirical = data,
    Lower_CI = lower_ci, Upper_CI = upper_ci
  ))
}




#-------------------------------------------------------------------------------
#' Quantile Function of the Generalized Gumbel Distribution (GGD)
#'
#' @param f Vector of non-exceedance probabilities (0 <= f <= 1).
#' @param para A parameter object returned by lh.parggd.
#' @return Quantiles for the specified probabilities.
#' @noRd
#'
quaggd <- function(f, para) {
  if (para$type != "ggd") {
    warning("Parameters are not of type 'ggd'.")
  }

  mu <- para$para[1]
  sigma <- para$para[2]
  h <- para$para[3]

  kap_para <- list(
    type = "kap",
    para = c(xi = mu, alpha = sigma, k = 0, h = h)
  )

  # ส่งต่อให้ lmomco คำนวณ
  return(lmomco::quakap(f, kap_para))
}





#-------------------------------------------------------------------------------
#' Quantile Function of the Three-Parameter Kappa Distribution (Fixed h)
#'
#' @param f Vector of non-exceedance probabilities (0 <= f <= 1).
#' @param para A parameter object returned by lh.park3d.hfix.
#' @return Quantiles for the specified probabilities.
#' @noRd
#'
quak3d.hfix <- function(f, para) {
  mu <- para$para[1]
  sigma <- para$para[2]
  k <- para$para[3]
  hfix <- para$para[4]

  kap_para <- list(
    type = "kap",
    para = c(xi = mu, alpha = sigma, k = k, h = hfix)
  )

  return(lmomco::quakap(f, kap_para))
}




#-------------------------------------------------------------------------------
#' Quantile Function of the Three-Parameter Kappa Distribution (Fixed k)
#'
#' @param f Vector of non-exceedance probabilities (0 <= f <= 1).
#' @param para A parameter object returned by lh.park3d.kfix.
#' @return Quantiles for the specified probabilities.
#' @noRd
#'
quak3d.kfix <- function(f, para) {
  mu <- para$para[1]
  sigma <- para$para[2]
  kfix <- para$para[3]
  h <- para$para[4]

  kap_para <- list(
    type = "kap",
    para = c(xi = mu, alpha = sigma, k = kfix, h = h)
  )

  return(lmomco::quakap(f, kap_para))
}




