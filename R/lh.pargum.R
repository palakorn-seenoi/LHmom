#------------------------------------------------------------------
# Estimating LH-me for Gumbel
# lh.pargum(data, eta=2)
#------------------------------------------------------------------

#' Estimate Parameters of the Gumbel Distribution using LH-moments
#'
#' This function estimates the parameters of the Gumbel distribution based on
#' the sample LH-moments. The Gumbel distribution is treated as a special case
#' of the Generalized Extreme Value (GEV) distribution with the shape parameter
#' (k) fixed at 0. The analytical estimation utilizes the first two LH-moments
#' and the Euler-Mascheroni constant.
#'
#' @param data A numeric vector of data values.
#' @param eta A non-negative integer (between 0 and 4) representing the order
#'   of the LH-moments. Default is 1.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{model}: The specific model name (\code{"gum"}).
#'   \item \code{para}: A named numeric vector containing the estimated parameters
#'     (\code{xi} for location, \code{alpha} for scale, and \code{k} = 0 for shape).
#'   \item \code{type}: The overarching distribution family (\code{"gev"}).
#'   \item \code{eta}: The order of the LH-moments used.
#'   \item \code{source}: The name of the function (\code{"lh.pargum"}).
#' }
#'
#' @export
lh.pargum = function(data,eta=1){

  L= lhmoms(data, eta=eta)$lambdas[1:2]

  sigma= 2*L[2]/( (eta+2)*( log(eta+2)-log(eta+1) ) )
  mu=    L[1]-sigma*(0.5772 + log(eta+1))

  para=c(mu,sigma,0)
  names(para) <- c("xi","alpha","k")

  return( list(model="gum", para=para, type="gev",
               eta=eta, source="lh.pargum"
  ) )
}
