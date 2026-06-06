
# for(eta in 0:4){
# #  cat("eta,para glo.ol=",eta,lh.parglo.k3h(data,eta)$para,"\n")
#   cat("eta,para glo.mk=",eta,lh.parglo(data,eta)$para,"\n","\n")
# #  cat("eta,para gpa.ol=",eta,lh.pargpa.k3h(data,eta)$para,"\n")
#   cat("eta,para gpa.mk=",eta,lh.pargpa(data,eta)$para,"\n","\n")
#   }
#------------------------------------------------
lh.parglo = function(data,eta=1){
  
  if(length(eta) != 1) stop("eta should be a scalar")
  if( eta %% 1 != 0 | eta < 0 | eta >= 5) {
    stop("eta must be a nonnegative integer less than 5")}
  
  samlh=lhmoms(data, eta=eta)
  t= samlh$ratios
  L=samlh$lambdas
  
  k= 3*((eta+2)^2)*t[3]-eta*(eta+3)
  k= -k/( (eta+4)*(eta+3) )
  
  num= 2*L[2]*k #* (-h)^(k+1)
  dem= (eta+2)* ( (eta+1)*beta(k+1,(eta+1) -k)
                  -(eta+2)*beta(k+1,(eta+2) -k) )
  sigma = num/dem
  fac= 1- ( (eta+1)*beta(k+1,(eta+1) -k) )
  mu= L[1]-fac*sigma/k

  para=c(mu, sigma, k)

  tau.diff= abs(lhmom.kap(c(para,-1),eta)$ratios[4]
                -t[4])
  
  t34.diff= ( abs(lhmom.kap(c(para,-1),eta)$ratios[3]
                  -t[3]) + tau.diff )
  
  return(list(model = "glo", para = para, eta=eta, 
              type="glo",tau.diff=tau.diff, ifail=0,
              t34.diff=t34.diff,source="lh.parglo"))
}
#------------------------------------------------
lh.pargpa= function(data,eta=1){
  

  if(length(eta) != 1) stop("eta should be a scalar")
  if( eta %% 1 != 0 | eta < 0 | eta >= 5) {
    stop("eta must be a nonnegative integer less than 5")}
  
  samlh=lhmoms(data, eta=eta)
  t= samlh$ratios
  L=samlh$lambdas
  
  k= (eta+3)*(1-3*t[3])/(3*t[3]+eta+3)
  
  num=2*L[2]*k 
  dem= (eta+2)*( (eta+1)*beta(k+1,(eta+1))
                 -(eta+2)*beta(k+1,(eta+2)) )
  sigma = num/dem
  fac= 1- ( (eta+1)*beta(k+1,(eta+1)) )
  mu= L[1]-fac*sigma/k
  
  para=c(mu, sigma, k)
  tau.diff= abs(lhmom.kap(c(para,1),eta)$ratios[4]
                -t[4])
  
  t34.diff= ( abs(lhmom.kap(c(para,1),eta)$ratios[3]
                  -t[3]) + tau.diff )
  
  return(list(model = "gpa", para = para, eta=eta, 
              type="gpa",tau.diff=tau.diff, ifail=0,
              t34.diff=t34.diff,source="lh.pargpa"))
}
# #------------------------------------------------
# lh.pargpa.k3h = function(data,eta=1){
#   
#   z=list()
#   z=lh.park3h(data,eta,hfix=1)
#   z$para= z$para[1:3]
#   z$type="gpa"
#   z$source="lh.pargpa"
#   z
# }
# #------------------------------------------------
# lh.parglo.k3h = function(data,eta=1){
#   
#   z=list()
#   z=lh.park3h(data,eta,hfix=-1)
#   z$para= z$para[1:3]
#   z$type="glo"
#   z$source="lh.parglo"
#   z
# }
