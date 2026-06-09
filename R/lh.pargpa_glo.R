
#------------------------------------------------------------------
# Estimating LH-me for glo, using Meshgi and Khalili (2009) SERRA
# lh.parglo(data, eta=2)
#------------------------------------------------------------------
#' LH-moments parameters estimation for GLO distribution
#'
#' Function to calculate LH-moments parameters estimation for GLO distribution
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

  return(list(para = para, eta=eta,
              type="glo", ifail=0, source="lh.parglo"))
}


#------------------------------------------------------------------
# Estimating LH-me for gpa, using Meshgi and Khalili (2009) SERRA
# lh.pargpa(data, eta=2)
#------------------------------------------------------------------
#' LH-moments parameters estimation for GPA distribution
#'
#' Function to calculate LH-moments parameters estimation for GPA distribution
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

  return(list(para = para, eta=eta,
              type="gpa",ifail=0,source="lh.pargpa"))
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
