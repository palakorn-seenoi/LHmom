
#------------------------------------------------------------------
# Estimating LH-me for glo, using Meshgi and Khalili (2009) SERRA
# lh.parglo(data, eta=2)
#------------------------------------------------------------------

#' Estimate Parameters of the Generalized Logistic (GLO) Distribution using LH-moments
#'
#' This function estimates the parameters of the Generalized Logistic (GLO)
#' distribution based on the sample LH-moments. The analytical formulas for
#' the parameter estimation are derived based on the methodology presented
#' by Meshgi and Khalili (2009).
#'
#' @param data A numeric vector of data values.
#' @param eta A non-negative integer (between 0 and 4) representing the order
#'   of the LH-moments. Default is 1.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{type}: The distribution type (\code{"glo"}).
#'   \item \code{para}: A named numeric vector containing the estimated parameters
#'     (\code{mu} for location, \code{sigma} for scale, \code{k} for shape).
#'   \item \code{eta}: The order of the LH-moments used.
#'   \item \code{source}: The name of the function (\code{"lh.parglo"}).
#'   \item \code{ifail}: A numeric indicator of success (0 for success).

#' }
#'
#' @references Meshgi, A., & Khalili, D. (2009). Comprehensive evaluation of
#' regional flood frequency analysis by L- and LH-moments. \emph{Stochastic
#' Environmental Research and Risk Assessment (SERRA)}, 23(1), 137-152.
#'
#' @export
lh.parglo = function(data,eta=1){

  if(length(eta) != 1) stop("eta should be a scalar")
  if( eta %% 1 != 0 | eta < 0 | eta >= 5) {
    stop("eta must be a nonnegative integer less than 5")}

  samlh=lhmoms(data, eta=eta)
  t= samlh$ratios
  L= samlh$lambdas

  k= 3*((eta+2)^2)*t[3]-eta*(eta+3)
  k= -k/( (eta+4)*(eta+3) )

  num= 2*L[2]*k #* (-h)^(k+1)
  dem= (eta+2)* ( (eta+1)*beta(k+1,(eta+1) -k)
                  -(eta+2)*beta(k+1,(eta+2) -k) )
  sigma = num/dem
  fac= 1- ( (eta+1)*beta(k+1,(eta+1) -k) )
  mu= L[1]-fac*sigma/k

  para=c(mu, sigma, k)
  names(para) <- c("mu","sigma","k")

  return(list(type="glo",
              para = para,
              eta=eta,
              source="lh.parglo",
              ifail=0))
}



#------------------------------------------------------------------
# Estimating LH-me for gpa, using Meshgi and Khalili (2009) SERRA
# lh.pargpa(data, eta=2)
#------------------------------------------------------------------

#' Estimate Parameters of the Generalized Pareto (GPA) Distribution using LH-moments
#'
#' This function estimates the parameters of the Generalized Pareto (GPA)
#' distribution based on the sample LH-moments. The analytical formulas for
#' the parameter estimation are derived based on the methodology presented
#' by Meshgi and Khalili (2009).
#'
#' @param data A numeric vector of data values.
#' @param eta A non-negative integer (between 0 and 4) representing the order
#'   of the LH-moments. Default is 1.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{type}: The distribution type (\code{"gpa"}).
#'   \item \code{para}: A named numeric vector containing the estimated parameters
#'     (\code{mu} for location, \code{sigma} for scale, \code{k} for shape).
#'   \item \code{eta}: The order of the LH-moments used.
#'   \item \code{source}: The name of the function (\code{"lh.pargpa"}).
#'   \item \code{ifail}: A numeric indicator of success (0 for success).

#' }
#'
#' @references Meshgi, A., & Khalili, D. (2009). Comprehensive evaluation of
#' regional flood frequency analysis by L- and LH-moments. \emph{Stochastic
#' Environmental Research and Risk Assessment (SERRA)}, 23(1), 137-152.
#'
#' @export
lh.pargpa= function(data,eta=1){

  if(length(eta) != 1) stop("eta should be a scalar")
  if( eta %% 1 != 0 | eta < 0 | eta >= 5) {
    stop("eta must be a nonnegative integer less than 5")}

  samlh=lhmoms(data, eta=eta)
  t= samlh$ratios
  L= samlh$lambdas

  k= (eta+3)*(1-3*t[3])/(3*t[3]+eta+3)

  num=2*L[2]*k
  dem= (eta+2)*( (eta+1)*beta(k+1,(eta+1))
                 -(eta+2)*beta(k+1,(eta+2)) )
  sigma = num/dem
  fac= 1- ( (eta+1)*beta(k+1,(eta+1)) )
  mu= L[1]-fac*sigma/k

  para=c(mu, sigma, k)
  names(para) <- c("mu","sigma","k")

  return(list(type="gpa",
              para = para,
              eta=eta,
              source="lh.pargpa",
              ifail=0))
}


# #------------------------------------------------
# lh.parglo = function(data,eta=1){
#
#   z=list()
#   z=lh.park3d.hfix(data,eta,hfix=-1)
#   z$para= z$para[1:3]
#   z$type="glo"
#   z$source="lh.parglo"
#   z
# }
# #------------------------------------------------
# lh.pargpa = function(data,eta=1){
#
#   z=list()
#   z=lh.park3d.hfix(data,eta,hfix=1)
#   z$para= z$para[1:3]
#   z$type="gpa"
#   z$source="lh.pargpa"
#   z
# }
