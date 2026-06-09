# R Code for LH-moments of GNO

# =============================================================
# Numerical Integration of Stable I_max function
# =============================================================

Imax <- function(n, theta) {

  integrand <- function(u) {
    (pnorm(u + theta)^(n - 1)) * dnorm(u)
  }

  trans= function(u){
    integrand(u/(1-u^2))*(1+u^2)/((1-u^2)^2)
  }
  # Integrate from -Inf to Inf
  # int_val <- integrate(integrand, lower = -Inf, upper = Inf)$value

  int_val = integrate(trans, lower = -1, upper = 1)$value
  # Integrate from -1 to +1

  # Multiply by the constant term from the stable formula
  return(n * exp((theta^2) / 2) * int_val)
}






# ==============================================================
# Function to calculate Theoretical LH-moments of GNO
# ==============================================================
# lhmom.gno(para=c(100,10,-0.3), eta=1)


#' Theoretical LH-moments function of GNO distribution
#'
#' Function to calculate theoretical LH-moments for GNO distribution
#' @export
lhmom.gno = function(para=NULL, eta=1){

  nmom=4
  if(length(eta) != 1) stop("eta should be a scalar")
  if( eta %% 1 != 0 | eta < 0 | eta >= 5) {
    stop("eta must be a nonnegative integer less than 5")}

  if(eta==0){
    z= lmomgno(vec2par(para,"gno"))
    z$eta=0
    return(z)
  }

  xi=para[1]
  alpha = para[2]
  k=para[3]

  lambdas= ratios= rep(NA, nmom)

  B1 = Imax(eta + 1, theta = -k)
  B2 = Imax(eta + 2, theta = -k)
  B3 = Imax(eta + 3, theta = -k)
  B4 = Imax(eta + 4, theta = -k)

  lambdas[1] <- xi + (alpha/k) * (1 - B1)

  #  lambdas[2] <- alpha*(eta+2)/(2*k) * (B1 -B2)
  lambdas[2] <- (eta+2)/2 * (-alpha/k) * (B2 - B1)

  #  lambdas[3] <- alpha*(eta+3)/(6*k) * ( -(eta+2)*B1
  #                                        + 2*(eta+3)*B2
  #                                        - (eta+4)*B3 )
  lambdas[3] <- (eta+3)/6 * (-alpha/k) * (     (eta+4)*B3
                                               - 2*(eta+3)*B2
                                               +   (eta+2)*B1  )

  #  lambdas[4] <- alpha*(eta+4)/(24*k) * ( (eta+2)*(eta+3)*B1
  #                                          - 3*(eta+3)*(eta+4)*B2
  #                                          + 3*(eta+4)*(eta+5)*B3
  #                                          -   (eta+5)*(eta+6)*B4 )
  lambdas[4] <- (eta+4)/24 * (-alpha/k) * (     (eta+6)*(eta+5)*B4
                                                - 3*(eta+5)*(eta+4)*B3
                                                + 3*(eta+4)*(eta+3)*B2
                                                -   (eta+2)*(eta+3)*B1  )



  ratios[2]=lambdas[2]/lambdas[1]
  for (j in 3:nmom) ratios[j]=lambdas[j]/lambdas[2]

  z=list()
  z$lambdas = lambdas;  z$ratios = ratios;  z$eta = eta
  z$type="gno"
  z
}




# =============================================================
# Sample LH-moments
# =============================================================

#' Sample LH-moments function
#'
#' Function to calculate sample LH-moments for sample data
#' @export
lhmoms = function(data, eta=NULL, nmom=5){

  if(is.null(eta) | max(eta) > 4 | min(eta) < 0) stop("eta should be given with positive and less than 5")
  if(all(eta %% 1 != 0) ) stop("eta should be an integer")

  z=list()
  n=length(data)
  x= sort(data)
  leta=length(eta)

  lambda= ratio= matrix(NA, nrow=nmom, ncol=5)

  for (ieta in eta){
    jeta=ieta+1

    for (r in 1:nmom){

      denom <- lchoose(n, ieta+r)
      ld = 0
      ld <- sum(sapply(seq(1, n), function(i) {
        wk <- sapply(seq(0, r - 1), function(k) {
          (-1)^k * exp(lchoose(r-1, k) + lchoose(i-1, ieta+r-1- k) + lchoose(n-i,k) - denom)
        })
        wk * x[i]
      }))

      lambda[r,jeta] <- ld/r

    } #end for r

    ratio[2,jeta]=lambda[2,jeta]/lambda[1,jeta]
    for (j in 3:nmom) ratio[j,jeta]=lambda[j,jeta]/lambda[2,jeta]

  } #end for ieta

  if(leta==1){
    z$eta=as.vector(eta); z$lambdas =as.vector(lambda[,jeta]);
    z$ratios = as.vector(ratio[,jeta])

  }else{
    z$eta=eta; z$lambdas = lambda; z$ratios = ratio
  }
  z
}



# =============================================================
# LHme for gno
# =============================================================
# lh.pargno(data,eta=4)

#' LH-moments estimation parameters for GNO distribution
#'
#' Function to calculate LH-moments estimation parameters for GNO distribution
#' @importFrom nleqslv nleqslv
#' @importFrom lmomco  pargno
#' @export
lh.pargno <- function(data, eta=NULL, opt=FALSE, ntry=5){

  z=list()
  if(length(eta) != 1) stop("eta should be a scalar")
  if( eta %% 1 != 0 | eta < 0 | eta >= 5) {
    stop("eta must be a nonnegative integer less than 5")}

  if(eta==0){
    z= pargno(lmoms(data))
    z$eta=0
    return(z)
  }

  samlh=lhmoms(data, eta=eta)
  t3= samlh$ratios[3]

  if(opt != TRUE){     # obtain LH-me without numerical optimization

    aa=matrix(NA,5,4)
    # aa[1,]= c(0,	-2.0462,	0,	-0.5016)
    aa[2,]= c(0.2825,	-2.5672,	-0.1521,	-0.6662)
    aa[3,]=	c(0.4720,	-2.9385,	-0.2351,	-0.9060)
    aa[4,]=	c(0.6142,	-3.2326,	-0.2694,	-1.1471)
    aa[5,]=	c(0.7276,	-3.4776,	-0.2746,	-1.3719)

    jeta=eta+1
    k= aa[jeta,1]+ aa[jeta,2]*t3 + aa[jeta,3]*t3^2 +aa[jeta,4]*t3^3

  }else if(opt==TRUE){

    init=initk(data, model="gno",ntry=ntry)

    itry=1
    mysol=list()

    #----------------------------------------------------
    obj.lhgno <- function(k, t3=t3, eta=eta) {

      if ( k < -1 ) return(10^6 )   # k for hosking style

      B1 = Imax(eta + 1, theta = -k)
      B2 = Imax(eta + 2, theta = -k)
      B3 = Imax(eta + 3, theta = -k)

      lambda3 <- (eta+3)/6 * ( - (eta+2)*B1 + 2*(eta+3)*B2
                               - (eta+4)*B3 )
      lambda3/((eta+2)/2 * (B1-B2) ) - t3
    }
    #----------------------------------------------------

    while(itry <= ntry){

      mysol = nleqslv(x= init[itry], fn=obj.lhgno, method="Broyden",
                      control=list(maxit=300, ftol=1e-6),
                      t3=t3, eta=eta)

      if(mysol$termcd <= 2 ){
        z$Ifail= 0
        k= mysol$x   # k for hosking style
        z$f= mysol$fvec
        break
      }else{
        itry=itry+1
        if(itry == ntry) {
          warning("no LHme solution found")
          z$Ifail= 5
          break }
      } #end if mysol
    } # end while

    # k <- tryCatch({
    #   uniroot(obj.fun.gno, interval = c(-0.99, 0.99))$root
    #   }, error = function(e) NA)
    #
    # if (is.na(k)) return(c(xi=NA, alpha=NA, k=NA))

  } # end if opt

  B1 = Imax(eta + 1, theta = -k)
  B2 = Imax(eta + 2, theta = -k)

  alpha <- 2*k*samlh$lambdas[2] / ( (eta+2) * (B1-B2) )

  xi <- samlh$lambdas[1] - (alpha / k) * (1-B1)

  para=c(xi, alpha, k)
  return(list(type = "gno", para = para, eta=eta, source="lh.pargno"))
}



#----------- External Packages ---------------------------------------

#library(lmomco)
# pargno function

#library(nleqslv)
# nleqslv function

#----------- sample running ---------------------------------------


#eta=0
#pargno(lmoms(data))
#lh.pargno(data,eta=eta, opt=FALSE)
#lh.pargno(data,eta=eta, opt=TRUE)
#
#for (eta in 0:4){
#  lhold= lh.pargno(data,eta=eta, opt=FALSE)
#  lhnew= lh.pargno(data,eta=eta, opt=TRUE)
#  cat("eta, para.old,para.new=",eta,lhold$para,lhnew$para,"\n")
#}

