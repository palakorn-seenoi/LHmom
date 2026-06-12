#---------------------------------------------
#' Generate Initial Parameters for LH-moments Optimization
#'
#' An internal helper function to generate a vector of initial estimates for
#' the shape parameter (typically \code{k}) used in the numerical optimization
#' of LH-moments. It uses standard L-moments to find a baseline parameter
#' and generates variations around it to prevent the solver from getting stuck.
#'
#' @param data A numeric vector of data values.
#' @param model A character string specifying the distribution model abbreviation
#'   (e.g., \code{"gev"}, \code{"gno"}). This string is used to dynamically call
#'   the corresponding \code{lmomco} function (e.g., \code{pargev}).
#' @param ntry An integer specifying the number of initial parameter combinations
#'   to generate. Minimum is 5. Default is 5.
#'
#' @return A numeric vector of length \code{max(5, ntry)} containing various
#'   initial values for the shape parameter.
#'
#' @importFrom lmomco  lmoms
#' @import stats
#' @keywords internal
#' @noRd
initk = function(data, model=NULL, ntry=5){

  init= rep(NA,max(5,ntry))
  lmom0=lmoms(data)

  parmodel = paste("par",model,sep="")
  init[1]= match.fun(parmodel)(lmom0)$para[3]
  init[2]= -0.25
  init[3]= 0.01
  init[4]= 0.25

  for (itry in 5:max(5,ntry)){
    init[itry]= init[1]+ (runif(1)-0.5)*0.3
  }

  init[which(init< -0.5)] = -0.49
  init[which(init > 0.5)] = 0.49
  init[which(abs(init) < 1e-2)] = -0.1

  init
}
