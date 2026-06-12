# for (eta in 0:4){
#   z=lh.parpe3(data,eta)
#   cat("eta, para=",eta,z$para,"\n")
# }

# -----------------------------------------------------------------
# Estimating LH-me for pe3
# lh.parpe3(data, eta=2)
#------------------------------------------------------------------
#' Estimate Parameters of the Pearson Type III (PE3) Distribution using LH-moments
#'
#' This function estimates the parameters of the Pearson Type III (PE3) distribution
#' based on the sample LH-moments. It provides two estimation methods: using
#' predefined polynomial coefficients (Log-Log Split Degree 5 based on sample LH-skewness)
#' or numerical optimization via \code{nleqslv}. If the numerical solver fails to
#' converge, the function falls back to the ordinary L-moments estimation (\code{eta = 0}).
#'
#' @param data A numeric vector of data values.
#' @param eta A non-negative integer (between 0 and 4) representing the order
#'   of the LH-moments. Default is 1.
#' @param opt A logical value. If \code{FALSE} (default), it estimates parameters
#'   using polynomial approximations. If \code{TRUE}, it utilizes numerical optimization.
#' @param ntry An integer specifying the maximum number of initialization attempts
#'   for the numerical optimization solver when \code{opt = TRUE}. Default is 10.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{para}: A named numeric vector containing the estimated parameters
#'     (\code{mu} for location, \code{sigma} for scale, \code{gamma} for shape).
#'   \item \code{eta}: The order of the LH-moments used.
#'   \item \code{type}: The distribution type (\code{"pe3"}).
#'   \item \code{ifail}: A numeric indicator of the optimization solver's status
#'     (0 for success, 2 for partial success, 5 for failure).
#'   \item \code{precision}: The final function value (\code{fvec}) from the solver.
#'   \item \code{source}: The name of the function (\code{"lh.parpe3"}).
#' }
#'
#' @importFrom lmomco  lmoms
#' @importFrom lmomco  parpe3
#' @importFrom nleqslv nleqslv
#' @export
lh.parpe3 = function(data, eta=1, opt=FALSE, ntry=10){

  z=list()
  if(length(eta) != 1) stop("eta should be a scalar")
  if( eta %% 1 != 0 | eta < 0 | eta >= 5) {
    stop("eta must be a nonnegative integer less than 5")}

  small=1e-5
  para=rep(NA,3)

  samlh=lhmoms(data, eta=eta)
  t3= samlh$ratios[3]

  if(eta==0){
    z= parpe3(lmoms(data))
    z$eta=0
    z$ifail=0
    z$source="lh.parpe3"
    return(z)
  }

  if(opt==FALSE){

    Lc <- matrix(
      c(-2.473889,   -2.752691 ,   -0.7396794 ,   -0.3278262, -6.890299e-02, -5.603958e-03,
        -41.971957,  -140.701142,  -194.8289081,  -136.8914578, -4.787208e+01, -6.724703e+00,
        -219.093042,  -821.415756, -1240.5409824,  -938.4834169, -3.544089e+02, -5.362215e+01,
        -579.087519, -2268.062122, -3563.6808333, -2801.0955085, -1.100065e+03, -1.730024e+02,
        -1104.623797, -4439.087182, -7148.1813554, -5756.9392236, -2.317581e+03, -3.735564e+02),
      nrow=5, ncol=6, byrow=TRUE)

    Hc = matrix(
      c(-3.967888,  -10.08536,  -16.75333,   -19.15619,   -11.59455,   -2.874049,
        14.457278,   93.82595,  176.05211,   106.10158,   -18.19470,  -27.208166,
        3.124592,  -74.63410, -492.98557, -1017.88881,  -876.52818, -273.316107,
        -4.821293, -157.37229, -750.22725, -1360.16706, -1081.11868, -318.053462,
        -2.823767, -137.68141, -654.05667, -1153.25357,  -888.00909, -253.368557),
      nrow=5, ncol=6, byrow=TRUE)

    at3 <- abs(t3)

    if (at3 < 1e-6) {
      alpha <- 10000
    } else if (at3 >= 1.0) {
      stop("Sample Tau3 (LH-skewness) from -1 to 1")
    } else {
      lt3 <- log(at3)
      cutoff <- 1/3

      vec.t3=c(1, lt3, lt3^2, lt3^3, lt3^4, lt3^5)

      # Log-Log Split Degree 5
      if (at3 <= cutoff) {
        ln_alpha <- Lc[eta+1,] %*% vec.t3
      } else {
        ln_alpha <- Hc[eta+1,] %*% vec.t3
      }
      alpha <- exp(ln_alpha)
    } # end if at3

    ifail= 0 # new add
    fvec = 1e-7 # new add-2

  }else if(opt==TRUE){

    init= abs(initk(data, model="pe3",ntry=ntry))
    init[1] = 4/(parpe3(lmoms(data))$para[3]^2)

    # ----- solve using nleqslv  -----
    itry = 1
    ifail = 10
    mysol = NULL

    ftol=ftolz= 1e-5
    xtol=xtolz= 1e-4

    #----------------------------------------------------
    obj.lhpe3 <- function(x, t3=t3, eta=eta) {

      tau3= abs(lhmom.pe3(c(0,x,1),eta)$ratios[3])   # alpha=x

      if( any(is.na(tau3)) ) return(10^6)
      return(tau3-abs(t3))
    }
    #----------------------------------------------------

    while (itry <= ntry) {

      mysol = nleqslv(x = init[itry],
                      fn = obj.lhpe3, method="Broyden",
                      control = list(maxit = 1000,
                                     ftol =ftol, xtol =xtol,
                                     allowSingular = TRUE),
                      t3=t3, eta=eta)

      if (abs(mysol$fvec) < ftol | mysol$termcd ==1){
        ifail = 0
        alpha = mysol$x
        fvec = mysol$fvec
        break

      }else{
        itry = itry + 1
        ftol= itry*ftolz
        xtol= itry*xtolz
        ifail=5  # new add
        fvec = 10^6  # new add-2

        if(mysol$termcd==2){
          ifail=2
          sol=mysol$x
          fvec=mysol$fvec
        }
        if (itry > ntry) {
          if(ifail != 2) ifail = 5
          break
        }
      } #end if sum
    } #end while

    if (ifail != 0) {
      if(ifail==2){
        alpha=sol
      }else{
        z$mysol = mysol
        cat("failure to solve LHme for pe3, at eta=",
            eta,"\n")
        cat("We use parpe3 at eta=0","\n")
        z=parpe3(lmoms(data))
        para=z$para
        names(para) <- c("mu", "sigma", "gamma")

        return(list(para = para, eta=eta,
                    type="pe3", ifail=ifail,
                    precision=fvec, source="lh.parpe3"))
      }
    } # end if

  } # end if opt

  L= samlh$lambdas[1:2]
  PWM=rep(NA,2)
  if(alpha <= 0) warning("alpha_hat should be positive")
  if(alpha < small) alpha= small

  for(r in 1:2){
    PWM[r]= cal_Lam(alpha, eta=eta, r=r)
  }

  be= L[2]/PWM[2]
  xi= L[1]- be*PWM[1]

  mu= xi+alpha*be
  sigma= be*sqrt(alpha)
  gam =2*sign(t3)/sqrt(alpha)

  # be  <- sigma / sqrt(alpha)
  # xi  <- mu - (alpha * be)

  # mu=L[1]
  # sigma= L[2]*sqrt(pi)*sqrt(alpha)*gamma(alpha)/gamma(alpha+0.5)
  # gam =2*sign(t3)/sqrt(alpha)

  para=c(mu, sigma, gam)
  names(para) <- c("mu", "sigma", "gamma")

  return(list(para = para, eta=eta,
              type="pe3", ifail=ifail,
              precision=fvec, source="lh.parpe3"))
}



# ==============================================================================
# Function: lhmom.pe3
# Description: Calculates Theoretical LH-moments (L1-L4) and LH-moment ratios (Tau2-Tau4)
#              for the Pearson Type III (PE3) distribution.
# Arguments:
#   para: A vector of length 3 c(xi, alpha, beta) OR a 'vec2par' object from lmomco.
#   eta:  The level of LH-moments (eta = 0 corresponds to standard L-moments).
# ==============================================================================

#' Calculate Theoretical LH-moments for the Pearson Type III (PE3) Distribution
#'
#' This function computes the theoretical LH-moments and LH-moment ratios for the
#' Pearson Type III (PE3) distribution. It supports parameter inputs either as a
#' \code{vec2par} object from the \code{lmomco} package or as a standard numeric vector.
#'
#' @param para A list containing a \code{vec2par} object (with \code{mu, sigma, gamma})
#'   OR a numeric vector of length 3 \code{c(xi, alpha, beta)}.
#' @param eta A non-negative integer representing the order of the LH-moments.
#'   Default is 0 (which corresponds to the standard L-moments).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{lambdas}: A named numeric vector of the first four theoretical LH-moments.
#'   \item \code{ratios}: A named numeric vector of the corresponding LH-moment ratios.
#'   \item \code{trim}: The trim level (fixed at 0 for LH-moments compatibility).
#'   \item \code{type}: The distribution type (\code{"pe3"}).
#'   \item \code{eta}: The order of the LH-moments used.
#'   \item \code{source}: The name of the function (\code{"lhmom.pe3"}).
#' }
#'
#' @importFrom lmomco vec2par
#' @export
lhmom.pe3 <- function(para, eta = 0) {

  # --- 1. Parameter Handling and Validation ---
  if (is.list(para) && !is.null(para$type) && para$type == "pe3") {
    # Case 1: Input is a 'vec2par' object from the lmomco package.
    # Convert (mu, sigma, gamma) to the (xi, alpha, beta) parameterization.

    mu    <- para$para[1]
    sigma <- para$para[2]
    gam <- para$para[3]

    alpha <- 4 / (gam^2)
    be  <- sigma / sqrt(alpha)
    xi    <- mu - (alpha * be)

  } else if (is.vector(para) && length(para) == 3) {
    # Case 2: Input is a standard numeric vector c(xi, alpha, beta).

    xi    <- para[1]
    alpha <- para[2]
    be  <- para[3]

  } else {
    stop("Invalid input: 'para' must be a vec2par object (pe3) or a vector c(xi, alpha, beta)")
  }

  if (eta < 0) stop("eta must be non-negative")

  # --- 2. Function to Calculate Generalized PWM (Beta_S) ---
  cal_beta <- function(S, alpha, be, xi) {

    if (abs(S) < 1e-9) return(xi + (alpha * be))
    term1 <- xi / (S + 1)

    res = cal_Omega(alpha,S)

    term2 <- be * res
    return(term1 + term2)
  }

  # --- 3. Calculate Lambda_1 to Lambda_4 ---
  L <- numeric(4)     # Initialize an empty numeric vector of length 4

  for (r in 1:4) {
    lambda_sum <- 0

    for (k in 0:(r - 1)) {
      # Combinatorial coefficients based on generalized Shifted Legendre polynomials
      term_factor <- (eta + r) / r
      term_sign   <- (-1)^(r - 1 - k)
      term_comb1  <- choose(r - 1, k)
      term_comb2  <- choose((r - 1) + k + eta, k + eta)

      coeff <- term_factor * term_sign * term_comb1 * term_comb2

      # Call the function to compute Beta_S (where S = eta + k)
      beta_val <- cal_beta(S = eta + k, alpha, be, xi)

      # cal_Omega(alpha,S = eta + k) =
      # cal_beta(S = eta + k, alpha, be=1, xi=0)

      # Accumulate the weighted sum
      lambda_sum <- lambda_sum + (coeff * beta_val)
    }

    L[r] <- lambda_sum
  }

  # --- 4. Calculate LH-moment Ratios (Tau2, Tau3, Tau4) ---
  Tau2 <- L[2] / L[1]  # L-CV (Coefficient of L-variation)
  Tau3 <- L[3] / L[2]  # L-Skewness
  Tau4 <- L[4] / L[2]  # L-Kurtosis

  # --- 5. Format Output (Matching 'lmomco' package structure) ---
  z <- list()

  z$lambdas <- L
  names(z$lambdas) <- c("L1", "L2", "L3", "L4")

  # The first element is conventionally NA since a ratio of L1/L1 is trivial
  z$ratios <- c(NA, Tau2, Tau3, Tau4)
  names(z$ratios) <- c("NA", "Tau2", "Tau3", "Tau4")

  z$trim <- 0
  z$type="pe3"
  z$eta <- eta
  z$source <- "lhmom.pe3"

  return(z)
}




#-------------------------------------------------
# Omega function for PE3 distribution
#-------------------------------------------------
#' Calculate the Omega Function for the PE3 Distribution
#'
#' An internal helper function to calculate the generalized Probability Weighted
#' Moments (PWM) component, Omega, for the Pearson Type III (PE3) distribution
#' via numerical integration.
#'
#' @param alpha A numeric shape parameter for the Gamma distribution component.
#' @param S A numeric value representing the moment order step (\code{eta + k}). Default is 1.
#'
#' @return A numeric value resulting from the numerical integration.
#' @keywords internal
cal_Omega <- function(alpha=NULL, S=1) {

  integrand <- function(t) {
    val <- numeric(length(t))
    val[t == 0] <- 0
    idx <- t > 0
    if (any(idx)) {
      log_val <- log(t[idx]) +
        dgamma(t[idx], shape=alpha, rate=1, log=TRUE) +
        (S * pgamma(t[idx], shape=alpha, rate=1, log.p=TRUE))
      val[idx] <- exp(log_val)
    }
    return(val)
  }

  # Integration Bounds
  lower_b <- qgamma(1e-8, shape = alpha, rate = 1)
  upper_b <- qgamma(1 - 1e-8, shape = alpha, rate = 1)

  if(lower_b < 0) lower_b <- 0

  res <- tryCatch({
    integrate(integrand, lower=lower_b, upper=upper_b,
              subdivisions = 2000,
              rel.tol = 1e-10, abs.tol = 1e-10)$value
  }, error = function(e) {
    return(NA)
  })

  return(res)
}



#----------------------------------------------------------
# Lamdbda function for PE3 distribution
#----------------------------------------------------------
#' Calculate the Lambda Function for the PE3 Distribution
#'
#' An internal helper function to compute the weighted sum of the Omega function
#' values, calculating the specific Lambda components for the PE3 parameter estimation.
#'
#' @param alpha A numeric shape parameter.
#' @param eta A non-negative integer representing the order of the LH-moments. Default is 1.
#' @param r An integer representing the moment order index. Default is 1.
#'
#' @return A numeric value representing the Lambda component.
#' @keywords internal
cal_Lam <- function(alpha=NULL, eta=1, r=1){

  Lam_val <- 0

  for (k in 0:(r - 1)) {
    term_factor <- (eta + r) / r
    term_sign   <- (-1)^(r - 1 - k)
    term_comb1  <- choose(r - 1, k)
    term_comb2  <- choose((r-1) + k + eta, k + eta)

    coeff <- term_factor * term_sign * term_comb1 * term_comb2

    Omega_val <- cal_Omega(alpha, S = eta + k)

    Lam_val <- Lam_val + coeff * Omega_val
  }
  return(Lam_val)
}
