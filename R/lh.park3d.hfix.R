
#for (eta in 0:4){
#  cat("eta,par=",eta,lh.park3d.hfix(data,eta, hfix=-.5)$para,"\n")
#}

# --------------------------------------------------------
# Estimate LHme for K3d with hfix
# lh.park3d.hfix(data,eta=2,hfix=-0.5)
#---------------------------------------------------------

#' Estimate Parameters of the Three-Parameter Kappa Distribution with Fixed h
#'
#' This function estimates the parameters of the three-parameter Kappa distribution
#' based on the sample LH-moments, given a fixed value for the shape parameter \code{h}.
#' The function utilizes numerical optimization (\code{nleqslv}) to solve for the
#' shape parameter \code{k}. If the numerical solver fails to converge, the function
#' implements a fallback mechanism by adopting the \code{k} parameter estimated from
#' the Generalized Extreme Value (GEV) distribution.
#'
#' @param data A numeric vector of data values.
#' @param eta A non-negative integer (between 0 and 4) representing the order
#'   of the LH-moments. Default is 1.
#' @param hfix A numeric scalar representing the fixed value for the shape parameter
#'   \code{h}. Default is 0.
#' @param ntry An integer specifying the maximum number of initialization
#'   attempts for the numerical optimization solver. Default is 10.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{para}: A named numeric vector containing the estimated parameters
#'     (\code{mu} for location, \code{sigma} for scale, \code{k} for shape, and \code{hfix}).
#'   \item \code{eta}: The order of the LH-moments used.
#'   \item \code{hfix}: The fixed value of the shape parameter \code{h} used in the estimation.
#'   \item \code{type}: The distribution type (\code{"kap"}).
#'   \item \code{ifail}: A numeric indicator of the optimization solver's status
#'     (0 for success, 2 for partial success/limit reached, 5 for failure).
#'   \item \code{precision}: The final function value (\code{fvec}) from the solver,
#'     indicating the precision of the root found.
#'   \item \code{source}: The name of the function (\code{"lh.park3d.hfix"}).
#' }
#'
#' @importFrom nleqslv nleqslv
#' @export
lh.park3d.hfix = function(data, eta=1, hfix=0,
                          ntry=10){

  z=list()
  if(length(eta) != 1) stop("eta should be a scalar")
  if( eta %% 1 != 0 | eta < 0 | eta >= 5) {
    stop("eta must be a nonnegative integer less than 5")}

  h=hfix
  if( h < -eta -1 | h > 5) stop("hfix out of bound")
  small= 1e-6

  if(h==0) hfix= h= -small
  if(h != 0 & abs(h) < small) hfix=h=sign(h)*small

  samlh= lhmoms(data, eta=eta)
  t3= samlh$ratios[3]

  # if(eta==0){
  #   z= park3d.hfix(lmoms(data), hfix=hfix)
  #   z$eta=0
  #   z$hfix=hfix
  #
  #   z$source="lh.park3d.hfix"
  #   return(z)
  # }

  init= initk(data, model="gev",ntry=ntry)

  #----------------------------------------------------
  obj.lhk3h <- function(x, t3=t3, eta=eta, hfix=hfix) {

    #  k=x
    if ( x <= -1 ) return(10^6 )
    if ( hfix < 0 & (x*hfix) <= -eta-1) return(10^6)
    #  if( k+0.725*h <= -1) return(c(10^6,10^6) )

    tau=rep(NA,4)
    tau[3:4]= lhmom.kap(c(0,1,x,hfix),eta)$ratios[3:4]

    if( any(is.na(tau[3:4])) ) return(10^6)
    if( tau[4] < (5*tau[3]^2-1)/4 ) return(10^6)

    return(tau[3]-t3)
  }
  #----------------------------------------------------
  # ----- solve using nleqslv  -----
  itry = 1
  ifail = 10
  mysol = NULL
  para=rep(NA,4)
  ftol=ftolz= small
  xtol=xtolz= 10*small

  while (itry <= ntry) {

    mysol = nleqslv(x = init[itry],
                    fn = obj.lhk3h, method="Broyden",
                    control = list(maxit = 1000,
                                   ftol =ftol, xtol =xtol,
                                   allowSingular = TRUE),
                    t3=t3, eta=eta, hfix=hfix)

    if (abs(mysol$fvec) < ftol| mysol$termcd ==1){
      ifail = 0
      para[3] = mysol$x
      fvec = mysol$fvec
      break

    }else{
      itry = itry + 1
      ftol=itry*ftolz
      xtol=itry*xtolz

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
      para[3]=sol
    }else{
      z$mysol = mysol
      cat("failure to solve LHme for K3D, at eta=",
          eta,"\n")
      cat("We use khat of GEV, with h=hfix","\n")
      para[3]=lh.pargev(data,eta)$para[3]
    }
  }

  L= samlh$lambdas[1:2]
  k=para[3]
  para[4]=hfix

  if(h < -eta -1) h= -eta-1 + small
  if(k < -1) k=-0.99999

  if(para[4] <0){
    num=2*L[2]*k* (-h)^(k+1)
    dem= (eta+2)* ( (eta+1)*beta(k+1,(-eta-1)/h -k)
                    -(eta+2)*beta(k+1,(-eta-2)/h -k) )

    sigma = num/dem
    fac= 1- ( (eta+1)*beta(k+1,(-eta-1)/h -k) / ((-h)^(k+1)) )
    mu=  L[1]-fac*sigma/k

  }else if(para[4] >0){
    num= 2*L[2]*k* (h)^(k+1)
    dem= (eta+2)*( (eta+1)*beta(k+1,(eta+1)/h)
                   -(eta+2)*beta(k+1,(eta+2)/h) )
    sigma = num/dem
    fac= 1- ( (eta+1)*beta(k+1,(eta+1)/h) / (h)^(k+1) )
    mu=  L[1]-fac*sigma/k
  } # end if

  para=c(mu, sigma, k, h)
  names(para)= c("mu", "sigma", "k", "hfix")

  return(list(para = para, eta=eta,
              hfix=hfix, type="kap",
              ifail=ifail, precision=fvec,
              source="lh.park3d.hfix"))
}
