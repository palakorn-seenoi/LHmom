
#------------------------------------------------------------
# Calculate the theoretical LH moments for GLO, GPA, GGD,
# K3D.hfix, and k3d.kfix. It are based on those of K4D.
# -----------------------------------------------------------
# para=c(10,1,-.3)
# lhmom.gum(para,eta=0); lhmom.gum(para,eta=4)
# lhmom.glo(para,eta=0); lhmom.glo(para,eta=4)
# lhmom.gpa(para,eta=0); lhmom.gpa(para,eta=4)
# lhmom.ggd(para,eta=0); lhmom.ggd(para,eta=4)
# lhmom.k3d.hfix(para,eta=0,hfix=-.5); lhmom.k3d.hfix(para,eta=4,hfix=-.5)
# lhmom.k3d.hfix(para,eta=0); lhmom.k3d.hfix(para,eta=4)
# lhmom.k3d.kfix(para,eta=0,kfix=-.3); lhmom.k3d.kfix(para,eta=4,kfix=-.3)
# lhmom.k3d.kfix(para,eta=0); lhmom.k3d.kfix(para,eta=4)


#------------------------------------------------------------
#' Theoretical LH-moments for the Generalized Logistic (GLO) Distribution
#'
#' This function computes the theoretical LH-moments for the Generalized
#' Logistic (GLO) distribution. It is evaluated as a special case of the
#' four-parameter Kappa distribution with the shape parameter h set to -1.
#'
#' @param para A numeric vector of three parameters: c(mu, sigma, k).
#' @param eta The order of LH-moments (default is 1).
#' @return A list containing the calculated LH-moments and the distribution type ("glo").
#' @export
lhmom.glo = function(para=NULL, eta=1){

  # para = c( mu, sigma, k)
  z=list()
  z= lhmom.kap(para=c(para,-1), eta=eta)
  z$type="glo"
  return(z)
}


#------------------------------------------------------------
#' Theoretical LH-moments for the Generalized Pareto (GPA) Distribution
#'
#' This function computes the theoretical LH-moments for the Generalized
#' Pareto (GPA) distribution. It is evaluated as a special case of the
#' four-parameter Kappa distribution with the shape parameter h set to 1.
#'
#' @param para A numeric vector of three parameters: c(mu, sigma, k).
#' @param eta The order of LH-moments (default is 1).
#' @return A list containing the calculated LH-moments and the distribution type ("gpa").
#' @export
lhmom.gpa = function(para=NULL, eta=1){

  # para = c( mu, sigma, k)
  z=list()
  z= lhmom.kap(para=c(para,1), eta=eta)
  z$type="gpa"
  return(z)
}

#------------------------------------------------------------
#------------------------------------------------------------
#' Theoretical LH-moments for the Generalized Gumbel (GGD) Distribution
#'
#' This function computes the theoretical LH-moments for the Generalized
#' Gumbel (GGD) distribution. It is evaluated as a special case of the
#' four-parameter Kappa distribution with the shape parameter k set to 0.
#'
#' @param para A numeric vector of three parameters: c(mu, sigma, h).
#' @param eta The order of LH-moments (default is 1).
#' @return A list containing the calculated LH-moments and the distribution type ("ggd").
#' @export
lhmom.ggd = function(para=NULL, eta=1){

  # para = c( mu, sigma, h)
  z=list()
  z= lhmom.kap(para=c(para[1:2],0,para[3]), eta=eta)
  z$type="ggd"
  return(z)
}

#------------------------------------------------------------
#' Theoretical LH-moments for the Three-Parameter Kappa Distribution with Fixed h
#'
#' This function computes the theoretical LH-moments for the three-parameter
#' Kappa distribution where the shape parameter h is fixed.
#'
#' @param para A numeric vector of three parameters: c(mu, sigma, k).
#' @param eta The order of LH-moments (default is 1).
#' @param hfix The fixed numeric value for the shape parameter h (default is 0).
#' @return A list containing the calculated LH-moments and the distribution type ("kap").
#' @export
lhmom.k3d.hfix = function(para=NULL, eta=1, hfix=0){

  # para = c( mu, sigma, k)
  z=list()
  z= lhmom.kap(para=c(para,hfix), eta=eta)
  z$type="kap"
  return(z)
}

#------------------------------------------------------------
#' Theoretical LH-moments for the Three-Parameter Kappa Distribution with Fixed k
#'
#' This function computes the theoretical LH-moments for the three-parameter
#' Kappa distribution where the shape parameter k is fixed.
#'
#' @param para A numeric vector of three parameters: c(mu, sigma, h).
#' @param eta The order of LH-moments (default is 1).
#' @param kfix The fixed numeric value for the shape parameter k (default is 0).
#' @return A list containing the calculated LH-moments and the distribution type ("kap").
#' @export
lhmom.k3d.kfix = function(para=NULL, eta=1, kfix=0){

  # para = c( mu, sigma, h)
  z=list()
  z= lhmom.kap(para=c(para[1:2],kfix,para[3]), eta=eta)
  z$type="kap"
  return(z)
}

#------------------------------------------------------------
#' Theoretical LH-moments for Gumbel Distribution
#'
#' This function computes the theoretical LH-moments for the Gumbel distribution
#' by evaluating it as a special case of the Generalized Extreme Value (GEV)
#' distribution with the shape parameter set to zero.
#'
#' @param para A vector of parameters c(mu, sigma)
#' @param eta The order of LH-moments (default is 1).
#' @return A list containing the LH-moments and distribution type.
#' @export
lhmom.gum = function(para=NULL, eta=1){

  # para = c(mu, sigma)
  z=list()
  z= lhmom.gev(para=c(para[1:2],0), eta=eta)
  z$type="gev"
  return(z)
}
#------------------------------------------------------------
