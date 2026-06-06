#--------------------------------------------------------------------
park3h=function (lmom, hfix=NULL, ntry=10,...) 
{
  
  h.fix=hfix
  if(is.null(h.fix) | h.fix < -5) {
    stop("h should not be less than -1")}
  
  small=1e-5
  if(hfix == 0) h.fix = -small
  if(hfix != 0 & abs(hfix) <= small) h.fix= sign(hfix)*small
  ifail=0
  
  para <- rep(NA, 4)
  names(para) <- c("xi", "alpha", "kappa", "h")
  
  T3 <- lmom$ratios[3]
  
  init=initk(data, model="gev", ntry=ntry)
  
  # init= (1 - 3 * T3)/(1 + T3)
  # if(init <= -1) init = -0.01
  nin = length(init)
  for (i in 1:nin){
   if(h.fix < 0 & (init[i]* h.fix) <= -1) {
     init[i]= -0.01
  }}
  
  itry=1
  mysol=list()
  
  while( itry <= ntry){
    
  mysol = nleqslv(x=init[itry], fn=solk3h, method="Broyden", 
                  control = list(maxit = 1000,
                                 ftol = 1e-5, xtol = 1e-4,
                                 allowSingular = TRUE),
                  h.fix=h.fix, T3=T3)
  
  if(mysol$termcd == 1){ #| mysol$termcd == 2){
    mysolG= mysol$x
    ifail=0
    break
  }else{
    itry=itry+1
    if(itry > ntry){
      ifail=5
      warning("No solution found in park3d: we use pargev")
      break
    }
  }
  } # end while
  
  if(ifail==0){
   para[3] <- G <- mysolG
   if(para[3] <= -1) {
     return("No feasible solution: k is smaller than -1")}
  }else if(ifail != 0){
    para[3]= G= pargev(lmoms(data))$para[3]
  }
  para[4] <- H <- h.fix 
  
  if (H > 0) {
    HH= exp(log(H)*(1+G))
    U1 <- exp(lgamma(1/H)+lgamma(1+G) - lgamma(1/H + 1 + G))/HH
    U2 <- exp(lgamma(2/H)+lgamma(1+G) - lgamma(2/H + 1 + G))*2/HH
  }else {
    HH = exp(log(-H)*(1+G))
    U1 <- exp(lgamma(-1/H - G)+lgamma(1+G) - lgamma(-1/H + 1))/HH
    U2 <- exp(lgamma(-2/H - G)+lgamma(1+G) - lgamma(-2/H + 1))*2/HH
  }
  
  alam2= U1-U2
  para[2]= lmom$lambdas[2]*G/alam2
  para[1]= lmom$lambdas[1]- para[2]*(1-U1)/G
  
  return(list(type = "kap", para = para, source = "park3h",
              ifail = ifail)) 
}
#---------------------- 
solk3h = function(x, h.fix=h.fix, T3=T3){
  
  if(x <= -1) return(10^6)
  if(h.fix < 0 &  (x * h.fix) <= -1 ) return(10^6)

  pk4= vec2par(vec=c(0,1, x, h.fix), type="kap")
  
  return(lmomkap(pk4)$ratios[3]-T3)
}
#---------------------