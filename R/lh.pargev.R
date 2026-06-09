# library(lmomco)
# library(nleqslv)

# for (eta in 0:4){
#   cat("eta,par=",eta,lh.pargev(data,eta)$para,"\n")
# }

#------------------------------------------------------------------
# Estimating LH-me for gevd, based on Wang's paper(1997)
# lh.pargev(data, eta=2)
#------------------------------------------------------------------
#' LH-moments parameters estimation for GEV distribution
#'
#' Function to calculate LH-moments parameters estimation for GEV distribution
#' @importFrom lmomco  lmoms
#' @importFrom lmomco  pargev
#' @importFrom nleqslv nleqslv
#' @export
lh.pargev = function(data, eta=1, opt=FALSE, ntry=5){

  if(length(eta) != 1) stop("eta should be a scalar")
  if( eta %% 1 != 0 | eta < 0 | eta >= 5) {
    stop("eta must be a nonnegative integer less than 5")}

  samlh = lhmoms(data, eta=eta)
  t3= samlh$ratios[3]

  if(eta==0){
    para= pargev(lmoms(data))$para
    return(list(type = "gev", para = para, eta=0,
                source="lh.pargev", ifail=0))
  }

  if(opt != TRUE){    # obtain LH-me without numerical optimization,
    # but use Wang's numbers

    cc=matrix(NA,5,4)
    cc[1,]=c(0.2849, -1.8213, 0.8140, -0.2835)
    cc[2,]=c(0.4823, -2.1494, 0.7269, -0.2103)
    cc[3,]=c(0.5914, -2.3351, 0.6442, -0.1616)
    cc[4,]=c(0.6618, -2.4548, 0.5733, -0.1273)
    cc[5,]=c(0.7113, -2.5383, 0.5142, -0.1027)

    jeta= eta+1
    k= cc[jeta,1]+ cc[jeta,2]*t3 + cc[jeta,3]*t3^2 +cc[jeta,4]*t3^3
    ifail=0

  }else if(opt==TRUE){  # obtain LH-me with numerical optimization

    init=initk(data, model="gev",ntry=ntry)

    itry=1
    mysol=list()

    #-----------------------------------------------
    obj.lhgev = function(x, t3=t3, eta=eta){

      if ( x < -1 ) return(10^6 )   # x is k for Hosking style

      num= ( -(eta+4)*(eta+3)^(-x) +2*(eta+3)*(eta+2)^(-x)
             -(eta+2)*(eta+1)^(-x) )*(eta+3)/6
      tau3 = num/(( (eta+1)^(-x) - (eta+2)^(-x) )*(eta+2)/2 )
      tau3-t3
    }
    #------------------------------------------------

    while(itry <= ntry){

      mysol = nleqslv(x=init[itry], fn=obj.lhgev, method="Broyden",
                      control=list(maxit=300, ftol=1e-6),
                      t3=t3, eta=eta)

      if(mysol$termcd <= 2 ){
        ifail= 0
        k= mysol$x   # k for hosking style
        #        f= mysol$fvec
        break
      }else{
        itry=itry+1
        if(itry == ntry) {
          warning("no LHme solution found")
          ifail= 5
          break }
      } #end if mysol
    } # end while

  } #end if opt

  deno=   gamma(1+k)*(eta+2)*((eta+1)^(-k)-(eta+2)^(-k))
  alpha= 2*k*samlh$lambdas[2]/deno
  xi =   samlh$lambdas[1]-(alpha/k)*(1- gamma(1+k)*(eta+1)^(-k))

  para= c(xi, alpha, k)
  names(para) <- c("xi","alpha","k")

  return(list(type = "gev", para = para, eta=eta,
              source="lh.pargev", ifail=ifail))
}


#-----------------------------------------------------------------
# Calculate population LH-moments for GEV, from Wang's paper (1997)
# lhmom.gev(para,eta=2)
#-----------------------------------------------------------------
#' Population LH-moments function of GEV distribution
#'
#' Function to calculate population LH-moments for GEV distribution
#' @export
lhmom.gev = function(para, eta=NULL){

  if(is.null(eta) | max(eta) >= 5 | min(eta) < 0) {
    stop("eta should be nonnegative and less than 5")}
  if(all(eta %% 1 != 0) ) stop("eta should be a nonnegative integer")

  nmom=4
  xi=    para[1]
  alpha= para[2]
  k =    para[3]

  small=1e-5
  if(k == 0) k= para[3]= -small
  if(k != 0 & abs(para[3]) <= small) k= para[3]= sign(k)*small

  lambdas= ratios= rep(NA, nmom)

  #-----------------------------------------------
  lh.Beta = function(feta, para){

    xi=    para[1]
    alpha= para[2]
    k =    para[3]

    small=1e-5
    if(abs(k) < small) k=sign(k)*small

    xi + (alpha/k)*(1- gamma(1+k)*(feta+1)^(-k))
  }
  #-----------------------------------------------

  B0 = lh.Beta(feta=eta+0, para)
  B1 = lh.Beta(feta=eta+1, para)
  B2 = lh.Beta(feta=eta+2, para)
  B3 = lh.Beta(feta=eta+3, para)

  lambdas[1]= B0
  lambdas[2]= (B1-B0)*(eta+2)/2
  lambdas[3]= ((eta+4)*B2-2*(eta+3)*B1 +(eta+2)*B0)*(eta+3)/6
  lambdas[4]= ( (eta+6)*(eta+5)*B3 -3*(eta+5)*(eta+4)*B2
                +3*(eta+4)*(eta+3)*B1 -(eta+3)*(eta+2)*B0 )*(eta+4)/24

  ratios[2]= lambdas[2]/lambdas[1]
  for (j in 3:nmom) ratios[j]= lambdas[j]/lambdas[2]

  z=list()
  z$eta=eta
  z$lambdas = as.vector(lambdas)
  z$ratios  = as.vector(ratios)
  names(z$lambdas) <- paste0("LHmom-",1:nmom)
  names(z$ratios)  <- paste0("LHtau-",1:nmom)
  z
}
