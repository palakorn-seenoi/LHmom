
# for(eta in 0:4){
#   cat("eta,para k3d=",eta,lh.park3d.kfix(data,eta,
#                                          kfix=-0.)$para,"\n")
#   cat("eta,para ggd=",eta,lh.parggd(data,eta)$para,"\n","\n")
# }

#----------------------------------------------------------
# Compute LHme for ggd
# lh.parggd(data,eta=2)
#----------------------------------------------------------

#' Estimate Parameters of the Generalized Gumbel (GGD) Distribution using LH-moments
#'
#' This function estimates the parameters of the Generalized Gumbel (GGD) distribution
#' based on the sample LH-moments. It evaluates the GGD as a special case of the
#' three-parameter Kappa distribution where the shape parameter \code{k} is fixed at 0.
#'
#' @param data A numeric vector of data values.
#' @param eta A non-negative integer (between 0 and 4) representing the order
#'   of the LH-moments. Default is 1.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{para}: A named numeric vector containing the estimated parameters
#'     (\code{mu} for location, \code{sigma} for scale, and \code{h} for shape).
#'   \item \code{type}: The distribution type (\code{"ggd"}).
#'   \item \code{source}: The name of the function (\code{"lh.parggd"}).
#'   \item \code{eta}: The order of the LH-moments used.
#' }
#'
#' @export
lh.parggd = function(data,eta=1){

  z=list()
  z=lh.park3d.kfix(data,eta,kfix=0)
  z$para= z$para[c(1:2,4)]
  names(z$para)= c("mu","sigma","h")
  z$type="ggd"
  z$source="lh.parggd"
  z
}


#---------------------------------------------------------
# Compute LHme for k3d with kfix
# lh.park3d.kfix(data,eta=2,kfix=-0.2)
#---------------------------------------------------------

#' Estimate Parameters of the Three-Parameter Kappa Distribution with Fixed k
#'
#' This function estimates the parameters of the three-parameter Kappa distribution
#' based on the sample LH-moments, given a fixed value for the shape parameter \code{k}.
#' It utilizes numerical optimization (\code{nleqslv}) to solve for the second shape
#' parameter \code{h}. If the numerical solver fails to converge, the function
#' implements a fallback mechanism by adopting the \code{h} parameter estimated from
#' the four-parameter Kappa distribution (\code{lh.parkap}).
#'
#' @param data A numeric vector of data values.
#' @param eta A non-negative integer (between 0 and 4) representing the order
#'   of the LH-moments. Default is 1.
#' @param kfix A numeric scalar representing the fixed value for the shape parameter
#'   \code{k}. Default is 0.
#' @param hlow A numeric scalar representing the lower bound for the shape parameter
#'   \code{h}. If \code{NULL}, it defaults to \code{-eta - 1}.
#' @param ntry An integer specifying the maximum number of initialization
#'   attempts for the numerical optimization solver. Default is 10.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{type}: The distribution type (\code{"kap"}).
#'   \item \code{para}: A named numeric vector containing the estimated parameters
#'     (\code{mu} for location, \code{sigma} for scale, \code{kfix}, and \code{h} for shape).
#'   \item \code{eta}: The order of the LH-moments used.
#'   \item \code{source}: The name of the function (\code{"lh.park3d.kfix"}).
#'   \item \code{ifail}: A numeric indicator of the optimization solver's status
#'     (0 for success, 2 for partial success/limit reached, 5 for failure).
#'   \item \code{kfix}: The fixed value of the shape parameter \code{k} used in the estimation.
#'   \item \code{precision}: The final function value (\code{fvec}) from the solver,
#'     indicating the precision of the root found.

#' }
#'
#' @importFrom nleqslv nleqslv
#' @export
lh.park3d.kfix = function(data, eta=1, kfix=0,
                          hlow=NULL,ntry=10){

  z=list()
  if(is.null(hlow)) hlow= -eta -1
  small=1e-5

  if(length(eta) != 1) stop("eta should be a scalar")
  if( eta %% 1 != 0 | eta < 0 | eta >= 5) {
    stop("eta must be a nonnegative integer less than 5")}

  k=kfix
  if( k < -1 | k > 1) stop("kfix out of bound")

  if(k==0) kfix= k= -small
  if(k != 0 & abs(k) < small) kfix=k=sign(k)*small

  samlh=lhmoms(data, eta=eta)
  t3= samlh$ratios[3]

  init= initk(data, model="gev",ntry=ntry)

  #----------------------------------------------------
  obj.lhk3k <- function(x, t3=t3, eta=eta,
                        kfix=kfix, hlow=hlow) {

    #  h=x
    if ( x < hlow | x >5) return(10^6 )
    if ( x < 0 & (x*kfix) <= -eta-1) return(10^6)
    #  if( k+0.725*h <= -1) return(c(10^6,10^6) )

    tau=rep(NA,4)
    tau[3:4]= lhmom.kap(c(0,1,kfix,x),eta)$ratios[3:4]

    if( any(is.na(tau[3:4])) ) return(10^6)
    if( tau[4] < (5*tau[3]^2-1)/4 ) return(10^6)

    return(tau[3]-t3)
  }
  #----------------------------------------------------
  # ----- solve using nleqslv  -----
  itry = 1
  mysol = NULL
  para=rep(NA,4)
  ftol=ftolz=1e-5
  xtol=xtolz=1e-4

  while (itry <= ntry) {

    mysol = nleqslv(x = init[itry],
                    fn = obj.lhk3k, method="Broyden",
                    control = list(maxit = 1000,
                                   ftol =ftol, xtol = xtol,
                                   allowSingular = TRUE),
                    t3=t3, eta=eta, kfix=kfix,
                    hlow=hlow)

    if (abs(mysol$fvec) <= ftol | mysol$termcd ==1){
      ifail = 0
      para[4] = mysol$x
      fvec = mysol$fvec
      break

    }else{
      itry = itry + 1
      ftol= itry*ftolz
      xtol= itry*xtolz

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
      para[4]=sol
    }else{

      cat("failure to solve LHme for lh.park3d.kfix at eta=",
          eta,"\n")
      cat("We use hhat from lh.parkap","\n")
      para[4]=lh.parkap(data,eta)$para[4]
    }
  }

  L= samlh$lambdas[1:2]
  h= para[4]

  if(h < hlow) h= hlow+small
  if(k < -1) k=-0.99999

  if(h <0){

    if(abs(k) <= small){
      den= (eta+2)*( digamma((-eta-2)/h)
                     -digamma((-eta-1)/h) )
      sigma= 2*L[2]/den
      mu=  L[1]-sigma*( digamma((-eta-1)/h)
                        -digamma(1) ) -sigma*log(-h)
    }else{
      num= 2*L[2]*k* (-h)^(k+1)
      den= (eta+2)* ( (eta+1)*beta(k+1,(-eta-1)/h -k)
                      -(eta+2)*beta(k+1,(-eta-2)/h -k) )

      sigma = num/den
      fac= 1- ( (eta+1)*beta(k+1,(-eta-1)/h -k) / ((-h)^(k+1)) )
      mu=  L[1]-fac*sigma/k
    }

  }else if(h >0){

    if(abs(k) <= small){
      den=(eta+2)*( digamma((eta+2)/h +1)
                    -digamma((eta+1)/h +1) )
      sigma= 2*L[2]/den
      mu=  L[1]-sigma*( digamma((eta+1)/h +1)
                        -digamma(1) ) - sigma*log(h)
    }else{
      num=2*L[2]*k* (h)^(k+1)
      den= (eta+2)*( (eta+1)*beta(k+1,(eta+1)/h)
                     -(eta+2)*beta(k+1,(eta+2)/h) )
      sigma = num/den
      fac= 1- ( (eta+1)*beta(k+1,(eta+1)/h) / (h)^(k+1) )
      mu=  L[1]-fac*sigma/k
    }
  } # end if

  para=c(mu, sigma, k, h)
  names(para)= c("mu", "sigma", "kfix", "h")

  return(list(type="kap",
              para = para,
              eta=eta,
              source="lh.park3d.kfix",
              ifail=ifail,
              kfix=kfix,
              precision=fvec) )
}
