
# #------------------------------------------------------
# lhmom.gpa = function(para=NULL,eta=1){
#   lhmom.kap(c(para,1),eta)
# }
# lhmom.glo = function(para=NULL,eta=1){
#   lhmom.kap(c(para,-1),eta)
# }
#---------------------------------------------------------

#' Estimate Parameters of the Four-Parameter Kappa Distribution using LH-moments
#'
#' This function estimates the parameters of the four-parameter Kappa (K4D) distribution
#' based on the sample LH-moments.
#' The estimation methodology follows Murshed et al. (2014).
#' It utilizes numerical optimization (\code{nleqslv})
#' to simultaneously solve for the two shape parameters, \code{k} and \code{h}.
#' If the optimization fails or becomes unstable, the function provides fallback
#' mechanisms utilizing fixed-parameter Kappa estimations. When \code{eta = 0},
#' it defaults to the ordinary L-moments estimation via \code{lmomco::parkap}.
#'
#' @param data A numeric vector of data values.
#' @param eta A non-negative integer (between 0 and 4) representing the order
#'   of the LH-moments. Default is 1.
#' @param snap.tau4 A logical value indicating whether to snap the sample L-kurtosis
#'   (\code{tau4}) downward if it slightly exceeds the theoretical upper bound.
#'   Default is \code{TRUE}.
#' @param nudge.tau4 A small numeric value used to adjust \code{tau4} downward
#'   if \code{snap.tau4 = TRUE}. Default is \code{1e-5}.
#' @param hlow A numeric scalar representing the lower bound for the shape parameter
#'   \code{h}. If \code{NULL}, it defaults to \code{-eta - 1}.
#' @param ntry An integer specifying the maximum number of initialization attempts
#'   for the numerical optimization solver. Default is 10.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{type}: The distribution type (\code{"kap"}).
#'   \item \code{para}: A named numeric vector containing the estimated parameters
#'     (\code{mu} for location, \code{sigma} for scale, \code{k} and \code{h} for shapes).
#'   \item \code{eta}: The order of the LH-moments used.
#'   \item \code{source}: The name of the function (\code{"lh.parkap"}).
#'   \item \code{ifail}: A numeric indicator of the optimization solver's status
#'     (0 for success, 1 for fallback success, 5 for failure).
#'   \item \code{precision}: The final function values from the solver.
#'   \item \code{ifailtext}: A descriptive message regarding the estimation success or failure.
#' }
#'
#' @references Murshed, S., Seo, Y.A., Park, J.S. (2014).
#' LH-moment estimation of a four parameter kappa distribution with hydrologic applications.
#' \emph{Stochastic Environmental Research and Risk Assessment}, 28, 253-262.
#'
#' @importFrom lmomco  lmoms
#' @importFrom lmomco  parkap
#' @importFrom nleqslv nleqslv
#' @export
lh.parkap = function(data, eta=1, snap.tau4= TRUE,
                     nudge.tau4=1e-5, hlow= NULL,
                     ntry=10){

  z=list()
  if(length(eta) != 1) stop("eta should be a scalar")
  if( eta %% 1 != 0 | eta < 0 | eta >= 5) {
    stop("eta must be a nonnegative integer less than 5")}

  small= 1e-5
  if(is.null(hlow)) hlow= -eta-1

  samlh=lhmoms(data, eta=eta)
  t3= samlh$ratios[3]
  t4= samlh$ratios[4]

  if(eta==0){
    z= parkap(lmoms(data))
    z$ifail=0

    if(any(is.na(z$para)) & snap.tau4==FALSE){
      z$message ="Failure to estim k4D at eta=0"
      z$ifail=5
      return(z)
    }

    if(any(is.na(z$para)) & snap.tau4==TRUE) {
      z= parkap(lmoms(data), snap.tau4=TRUE,
                nudge.tau4=nudge.tau4)
      z$ifail= -1
      if(any(is.na(z$para))) {
        z$message ="Failure to estim k4D at eta=0"
        z$ifail=5
        return(z)
      }
    }

    z$type="kap";  z$eta=0
    z$source="lh.park4d"
    # z$snap.tau4 = snap.tau4
    # z$nudge.tau4= nudge.tau4

    z$ifailtext=NULL
    if(z$ifail== -1){
      z$ifailtext="para estim under snap.tau4 & nudge.tau4"
    }else if(z$ifail==0){
      z$ifailtext="Successful parameter estimation."
    }
    return(z)
  } # end if eta=0

  init=initkh.k4d(data, ntry=ntry)

  #--------------------------------------------------------
  obj.lhk4d <- function(x, t3=t3,t4=t4, eta=eta,
                        hlow=hlow) {
    k=x[1]
    h=x[2]

    if ( k <= -1 ) return(rep(10^6,2) )     # k for hosking style
    if (h < 0 & (k*h) <= -eta-1) return(c(10^6,10^6) )

    if(h < hlow | h > 5 ) return(c(10^6,10^6) )

    #  if( k+0.725*h <= -1) return(c(10^6,10^6) )

    tau=rep(NA,4)
    tau[3:4]= lhmom.kap(c(0,1,k,h),eta)$ratios[3:4]

    if( any(is.na(tau[3:4])) ) return(rep(10^6,2))

    return(c(tau[3]-t3,tau[4]-t4))
  }
  #--------------------------------------------------------
  # ----- solve using nleqslv  -----
  itry = 1
  mysol=list()
  para=rep(NA,4)
  ftol=ftolz= small
  xtol=xtolz= 10*small

  if(snap.tau4==TRUE){ subt= subtz= nudge.tau4
  ifail=6
  t4.up=(5*(t3^2)+1)/6
  if( t4 >= t4.up) {
    t4= t4.up - subt }
  }else{ subt= subtz=0}

  tryCatch({
    while (itry <= ntry) {

      mysol =  tryCatch(nleqslv(x =as.vector(init[itry, 1:2]),
                                fn = obj.lhk4d,
                                control = list(maxit = 1000,
                                               ftol =ftol, xtol =xtol,
                                               allowSingular = TRUE),
                                t3=t3, t4=t4, eta=eta, hlow=hlow) )

      if (sum(abs(mysol$fvec)) < 2*ftol| mysol$termcd ==1){
        ifail = 0
        para[3:4] = mysol$x
        fvec = mysol$fvec
        break

      }else if(mysol$termcd==2){
        fvec =mysol$fvec
        sol =mysol$x

        if(abs(fvec[1]) <= ftol*10 & abs(fvec[2]) > ftol*10){
          z2=list()
          z2= lh.park3d.kfix(data, eta, kfix=sol[1])

          if(z2$ifail==0) {
            ifail=1
            para=z2$para
            break
          }else{
            itry = itry + 1
            ftol= itry*ftolz
            xtol= itry*xtolz

            if (itry > ntry) { ifail =5; break }
          }

        }else if(abs(fvec[1]) > ftol*10 & abs(fvec[2]) <= ftol*10){
          z2=list()
          z2= lh.park3d.hfix(data, eta,hfix=sol[2])

          if(z2$ifail==0) {
            ifail=1
            para=z2$para
            break
          }else{
            itry = itry + 1
            ftol= itry*ftolz
            xtol= itry*xtolz

            if (itry > ntry) { ifail =5; break }
          }

        }else if(abs(fvec[1]) <= ftol*10 & abs(fvec[2]) <= ftol*10){
          ifail=0
          para[3:4] = mysol$x
          break

        }else{

          itry = itry + 1
          ftol= itry*ftolz
          xtol= itry*xtolz

          if (itry > ntry) { ifail =5; break }
        }

      }else if(mysol$termcd >= 3){
        itry = itry + 1
        ftol= itry*ftolz
        xtol= itry*xtolz
        fvec= c(10^6, 10^6)

        if (itry > ntry) { ifail =5; break }

      } #end if sum
    } #end while
  }) # end trycatch

  if (ifail ==5) {
    # cat("failure to solve LHme for K4D, eta=",eta,"\n")

    z=list()
    z$para=rep(NA,4)
    z$ifail =ifail
    z$type="kap"

    z$eta=eta
    # z$snap.tau4 = snap.tau4
    # z$nudge.tau4=subt
    z$ifailtext="failue for k4d"
    return(z)
  }

  L= samlh$lambdas[1:2]
  k= para[3]
  h= para[4]

  if(h < hlow) h= hlow +small
  if(k < -1) k=-0.99999

  if(ifail==0){

    if(para[4] <0){
      num=2*L[2]*k* (-h)^(k+1)
      den= (eta+2)* ( (eta+1)*beta(k+1,(-eta-1)/h -k)
                      -(eta+2)*beta(k+1,(-eta-2)/h -k) )

      sigma = num/den
      fac= 1- ( (eta+1)*beta(k+1,(-eta-1)/h -k) / ((-h)^(k+1)) )
      mu=  L[1]-fac*sigma/k

    }else if(para[4] >0){
      num=2*L[2]*k* (h)^(k+1)
      den= (eta+2)*( (eta+1)*beta(k+1,(eta+1)/h)
                     -(eta+2)*beta(k+1,(eta+2)/h) )
      sigma = num/den
      fac= 1- ( (eta+1)*beta(k+1,(eta+1)/h) / (h)^(k+1) )
      mu=  L[1]-fac*sigma/k
    } # end if

    para= c(mu, sigma, k, h)
    names(para)= c("mu","sigma","k","h")
  } # end ifail

  if(ifail==5){
    ifailtext="failue for k4d"
  }else if(ifail==0){
    ifailtext="Successful parameter estimation."
  }else if(ifail==1){
    ifailtext=paste("x within xtol, but fcn value is > ftol.",
                    "para are obtained using lh.park3d.kfix or lh.park3d.hfix")
  }

  return(list(type = "kap",
              para = para,
              eta=eta,
              source="lh.parkap",
              ifail=ifail,
              precision=fvec, #snap.tau4=snap.tau4, nudge.tau4=subt,
              ifailtext=ifailtext))
}




#--------------------------------------------------------------
# Calculate the theoretical LH moments for K4D,
# based on Murshed et al.(2014) SERRA
#--------------------------------------------------------------

#' Calculate Theoretical LH-moments for the Four-Parameter Kappa Distribution
#'
#' This function computes the theoretical LH-moments and LH-moment ratios for
#' the four-parameter Kappa (K4D) distribution given its parameters. The computations
#' are analytically derived based on the formulas presented by Murshed et al. (2014).
#'
#' @param para A numeric vector of four parameters: c(xi, alpha, k, h), corresponding
#'   to the location, scale, and the two shape parameters of the K4D distribution.
#' @param eta A non-negative integer (between 0 and 4) representing the order
#'   of the LH-moments. Default is 1.
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{lambdas}: A named numeric vector of the first four theoretical LH-moments.
#'   \item \code{ratios}: A named numeric vector of the corresponding LH-moment ratios.
#'   \item \code{eta}: The order of the LH-moments used.
#'   \item \code{type}: The distribution type (\code{"kap"}).
#' }
#'
#' @references Murshed, S., Seo, Y.A., Park, J.S. (2014).
#' LH-moment estimation of a four parameter kappa distribution with hydrologic applications.
#' \emph{Stochastic Environmental Research and Risk Assessment}, 28, 253-262.
#'
#' @importFrom lmomco   lmoms
#' @importFrom lmomco   vec2par
#' @importFrom lmomco   lmomkap
#' @export
lhmom.kap = function(para=NULL, eta=1){

  nmom=4
  if(length(eta) != 1) stop("eta should be a scalar")
  if( eta %% 1 != 0 | eta < 0 | eta >= 5) {
    stop("eta must be a nonnegative integer less than 5")}

  z=list()
  if(eta==0){
    z= lmomkap(vec2par(para,"kap"))
    z$eta=0
    return(z)
  }

  xi= para[1]
  alpha = para[2]
  k= para[3]
  h= para[4]

  if(h == 0) {
    z= lhmom.gev(para[1:3],eta)
    z$eta = eta
    z$type="kap"
    return(z)
  }

  small=1e-5
  if(k != 0 & abs(k) < small) k=sign(k)*small
  if(k == 0) k= -small
  if(h != 0 & abs(h) < small) h=sign(h)*small

  lambdas= ratios= rep(NA, nmom)
  B=rep(NA,8)

  if(k <= -1) stop("k should be > -1 in lhmom.k4d")

  if(h < 0){
    if( (k*h) <= -eta-1){
      stop("k*h should be > -eta-1 in lhmom.k4d")

    }else{
      for(ieta in 0:7){

        num= (ieta+1)
        den= ((-h)^(k+1))

        bb= num*beta(k+1,(-ieta-1)/h -k)/den
        B[ieta+1]= xi + (alpha/k)*(1-bb)
      }
    } #end if k

  }else if(h >0){
    for(ieta in 0:7){
      num= (ieta+1)
      den= (h^(k+1))

      bb= num*beta(k+1,(ieta+1)/h)/den
      B[ieta+1]= xi + (alpha/k)*(1-bb)
    }
  } # end if h

  lam=tau=rep(NA,4)
  jeta=eta+1
  lam[1]=B[jeta]
  lam[2]=(eta+2)*(B[jeta+1]-B[jeta])/2
  lam[3]=(eta+3)*( (eta+4)*B[jeta+2]-2*(eta+3)*B[jeta+1]
                   +(eta+2)*B[jeta] )/6
  lam[4]=(eta+4)*( (eta+6)*(eta+5)*B[jeta+3]
                   -3*(eta+5)*(eta+4)*B[jeta+2]
                   +3*(eta+4)*(eta+3)*B[jeta+1]
                   -(eta+3)*(eta+2)*B[jeta] )/24

  tau[2]=lam[2]/lam[1]
  tau[3]=lam[3]/lam[2]
  tau[4]=lam[4]/lam[2]

  z$lambdas = lam
  z$ratios = tau
  names(z$lambdas) <- paste0("LHmom-",1:nmom)
  names(z$ratios)  <- paste0("LHtau-",1:nmom)
  z$eta = eta
  z$type="kap"
  z
}




#--------------------------------------------------------
# Initial parameters of K4D distribution
#--------------------------------------------------------

#' Generate Initial Parameters for Four-Parameter Kappa Optimization
#'
#' An internal helper function used to initialize the shape parameters (\code{k} and \code{h})
#' for the numerical optimization routines involved in estimating the four-parameter
#' Kappa distribution via LH-moments.
#'
#' @param data A numeric vector of data values.
#' @param ntry An integer specifying the number of initial parameter combinations
#'   to generate. Minimum is 5. Default is 5.
#'
#' @return A matrix containing \code{ntry} rows and 2 columns, representing
#'   various starting points for the \code{k} and \code{h} shape parameters.
#' @keywords internal
initkh.k4d = function(data, ntry=5){

  init= matrix(NA,max(5,ntry),2)
  lmom0=lmoms(data)

  init[1,]= parkap(lmom0, snap.tau4=TRUE,
                   nudge.tau4=1e-3)$para[3:4]

  if(any(is.na(init[1,1:2]))) {
    init[1,]= c(pargev(lmom0)$para[3],-0.05) }

  init[2,]= c(-0.25,-0.75)
  init[3,]= lh.park3d.hfix(data,eta=0, hfix=-1)$para[3:4]
  init[4,]= lh.park3d.hfix(data,eta=0, hfix=1)$para[3:4]
  init[5,]= c(0.25, 0.75)

  if(ntry >= 6){
    for (itry in 6:max(5,ntry)){
      init[itry,]= init[1,]+ rep((runif(1)-0.5)*0.3,2)
    }}
  init[which(init[,1]< -0.5),1] = -0.49
  init[which(init[,1] > 0.5),1] = 0.49
  init[which(abs(init[,1]) < 1e-2),1] = 0.01
  init[which(init[,2]< -2),2] = -2
  init[which(init[,2] > 1.5),2] = 1.5
  init[which(abs(init[,2]) < 1e-2),2] = 0.01

  init
}

