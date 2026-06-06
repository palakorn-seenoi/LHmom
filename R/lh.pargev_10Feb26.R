library(lmomco)
library(nleqslv)

# for (eta in 0:4){
#   cat("eta,par=",eta,lh.pargum(data,eta)$para,"\n")
# }
# lh.pargum(data,eta=0)
# pargum(lmoms(data))

#---------------------------------------------  
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

#------------------------------------------------------------
lh.pargum = function(data,eta=1){
  
  L=lhmoms(data, eta=eta)$lambdas[1:2]
  sigma= 2*L[2]/( (eta+2)*( log(eta+2)-log(eta+1) ) )
  mu= L[1]-sigma*(0.5772 + log(eta+1))
  
  tau4= lhmom.gev(para=c(mu,sigma,0),eta)$ratios[4]
  tau.diff=abs(tau4-lhmoms(data, eta=eta)$ratios[4])
  
  tau3= lhmom.gev(para=c(mu,sigma,0),eta)$ratios[3]
  t34.diff= ( abs(tau3-lhmoms(data, eta=eta)$ratios[3])
              + tau.diff )
  
  return( list(model="gum", para=c(mu,sigma,0), type="gev",
              eta=eta, source="lh.pargum",
              tau.diff=tau.diff, t34.diff=t34.diff) )
}
#------------------------------------------------------------------
# Estimating LH-me for gevd
# lh.pargev(data, eta=2)
#------------------------------------------------------------------
lh.pargev = function(data, eta=1, opt=FALSE, ntry=5){  
  
  z=list()
  if(length(eta) != 1) stop("eta should be a scalar")
  if( eta %% 1 != 0 | eta < 0 | eta >= 5) {
    stop("eta must be a nonnegative integer less than 5")}
  
  samlh=lhmoms(data, eta=eta)
  t3= samlh$ratios[3]
  t4= samlh$ratios[4]
  
  if(eta==0){
    z= pargev(lmoms(data))
    z$eta=0
    z$ifail=0
    z$tau.diff=abs(lhmom.gev(z$para,eta)$ratios[4]-t4)
    z$t34.diff= ( abs(lhmom.gev(z$para,eta)$ratios[3]-t3)
                  + z$tau.diff )
    z$source="lh.pargev"
    return(z)
  }
  
  if(opt != TRUE){    # obtain LH-me without numerical optimization
    
    cc=matrix(NA,5,4)
    cc[1,]=c(0.2849, -1.8213, 0.8140, -0.2835)
    cc[2,]=c(0.4823, -2.1494, 0.7269, -0.2103)
    cc[3,]=c(0.5914, -2.3351, 0.6442, -0.1616)
    cc[4,]=c(0.6618, -2.4548, 0.5733, -0.1273)
    cc[5,]=c(0.7113, -2.5383, 0.5142, -0.1027)
  
    jeta=eta+1
    k= cc[jeta,1]+ cc[jeta,2]*t3 + cc[jeta,3]*t3^2 +cc[jeta,4]*t3^3
    z$ifail=0
    
  }else if(opt==TRUE){  # obtain LH-me with numerical optimization
    
    init=initk(data, model="gev",ntry=ntry)
    
    itry=1
    mysol=list()

    #-----------------------------------------------
    obj.lhgev = function(x, t3=t3, eta=eta){
      
      if ( x < -1 ) return(10^6 )   # x is k for hosking style
      
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
        z$ifail= 0
        k= mysol$x   # k for hosking style
        z$f= mysol$fvec
        break
      }else{
        itry=itry+1
        if(itry == ntry) {
          warning("no LHme solution found")
          z$ifail= 5
          break }
      } #end if mysol
    } # end while
    
  } #end if opt

  dem= gamma(1+k)*(eta+2)*((eta+1)^(-k)-(eta+2)^(-k))
  alpha= 2*k*samlh$lambdas[2]/dem
  xi = samlh$lambdas[1]-(alpha/k)*(1- gamma(1+k)*(eta+1)^(-k))

  z$eta=eta
  z$para=c(xi, alpha, k)
  
  z$tau.diff=abs(lhmom.gev(z$para,eta)$ratios[4]-t4)
  z$t34.diff= ( abs(lhmom.gev(z$para,eta)$ratios[3]-t3)
                + z$tau.diff )
  
  z$type='gev'; z$source="lh.pargev"
  z
}

#-----------------------------------------------------------------
# population LH-moments for GEV
# lhmom.gev(para,eta=2)
#-----------------------------------------------------------------
lhmom.gev = function(para, eta=NULL){  
  
  if(length(eta) != 1) stop("eta should be a scalar")
  if( eta %% 1 != 0 | eta < 0 | eta >= 5) {
    stop("eta must be a nonnegative integer less than 5")}
  
  nmom=4
  xi= para[1]
  alpha=para[2]
  k = para[3]
  
  small=1e-5
  if(k == 0) k=para[3]= -small
  if(k != 0 & abs(para[3]) <= small) k=para[3]= sign(k)*small
  
  lambdas= ratios= rep(NA, nmom)

  #----------------------------------------
  lh.Beta = function(feta, para){
    xi=para[1]
    alpha=para[2]
    k = para[3]
    
    small=1e-5
    if(abs(k) < small) k=sign(k)*small
    
    gam=gamma(1+k)
    xi + (alpha/k)*(1- gam*(feta+1)^(-k))
  }
  #----------------------------------------
  
  B0 = lh.Beta(feta=eta, para)
  B1 = lh.Beta(feta=eta+1, para)
  B2 = lh.Beta(feta=eta+2, para)
  B3 = lh.Beta(feta=eta+3, para)

  lambdas[1]= B0
  lambdas[2]= (B1-B0)*(eta+2)/2
  lambdas[3]= ((eta+4)*B2-2*(eta+3)*B1 +(eta+2)*B0)*(eta+3)/6
  lambdas[4]= ( (eta+6)*(eta+5)*B3 -3*(eta+5)*(eta+4)*B2 
               +3*(eta+4)*(eta+3)*B1 -(eta+3)*(eta+2)*B0 )*(eta+4)/24

  ratios[2]=lambdas[2]/lambdas[1]
  for (j in 3:nmom) ratios[j]=lambdas[j]/lambdas[2]

  z=list()
    z$eta=eta; z$lambdas =as.vector(lambdas)
    z$ratios = as.vector(ratios)
  z
}

#------------------------------------------------------------------
# sample LH-moments from data
#  lhmoms(data,eta=seq(0,4))
#------------------------------------------------------------------
lhmoms= function(data, eta=NULL, nmom=5){      
  
  if(is.null(eta) | max(eta) >= 5 | min(eta) < 0) {
    stop("eta should be nonnegative and less than 5")}
  if(all(eta %% 1 != 0) ) stop("eta should be a nonnegative integer")
  
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
        (-1)^k * ( exp(lchoose(r-1, k) + lchoose(i-1, ieta+r-1- k) 
                     + lchoose(n-i,k) - denom) )
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
    z$eta=eta; z$lambdas = t(lambda); z$ratios = t(ratio)
  }
  z
}
#----------------------------------------------------------------
# data=haenam[,2]
# wang.test.lhgev(data, eta=3)
#----------------------------------------------------------------
wang.test.lhgev =function(data, eta=2){
  
  # tau = population L-moments
  # t = sample L-moment
  # n= sample size
  
  n=length(data)
  LH= lh.pargev(data,eta=eta)   # computing LHme for gevd
  tau= lhmom.gev(LH$para,eta)$ratios  
  
  t= lhmoms(data,eta)$ratios  # sample L-moments and L-ratios
  
  # lh.pargev: computing LHme for gevd
  # lhmom.gev: population L-moments and L-ratios for gevd
  # lhmoms: sample L-moments and L-ratios up to eta=4

  w=list()
  bb=cc=matrix(NA,5,5)
  bb[1,]=c(0.0745,  0.0555, 0.0067, -0.3090, 0.2240)
  bb[2,]=c(0.0579, -0.0328, 0.1524, -0.4102, 0.2672)
  bb[3,]=c(0.0488, -0.0527, 0.1620, -0.3856, 0.2566)
  bb[4,]=c(0.0380, -0.0309, 0.0354, -0.1233, 0.0878)
  bb[5,]=c(0.0241,  0.0024,-0.0813,  0.0733,-0.0210)
  
  cc[1,]=c(1.0100, -0.0282, -2.9336, 4.0801, -1.0874)
  cc[2,]=c(1.3403, -0.8291, -3.8777, 9.5371, -5.7866)
  cc[3,]=c(1.8800, -2.2233, -2.5825, 10.435, -7.3887)
  cc[4,]=c(2.6784, -4.8418,  3.5255, 2.3736, -3.2076)
  cc[5,]=c(3.7793, -8.3485,  11.517, -7.9095, 1.9459)
  
  jeta=eta+1
  B= bb[jeta,1] +bb[jeta,2]*tau[3] +bb[jeta,3]*tau[3]^2 
  B= B +bb[jeta,4]*tau[3]^3 +bb[jeta,5]*tau[3]^4
  
  C= cc[jeta,1] +cc[jeta,2]*tau[3] +cc[jeta,3]*tau[3]^2 
  C= C +cc[jeta,4]*tau[3]^3 +cc[jeta,5]*tau[3]^4
  
  cond.sigma_t4 =sqrt(B/n + C/(n^2))
  
  w$eta= eta
  w$z.test= (t[4]-tau[4])/cond.sigma_t4
  
  w$p.value= 2*pnorm( abs(w$z.test),lower.tail=FALSE)
  w$cond.sigma= cond.sigma_t4
  w
}