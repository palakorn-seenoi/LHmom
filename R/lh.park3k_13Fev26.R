library(lmomco)
library(nleqslv)
#---------------------------------------------------------
lh.park3k = function(data, eta=1, kfix=0, 
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
  t4= samlh$ratios[4]
  
  init= initk(data, model="gev",ntry=ntry)
  
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

      cat("failure to solve LHme for lh.park3k at eta=",
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
      dem= (eta+2)*( digamma((-eta-2)/h)
                    -digamma((-eta-1)/h) )
      sigma= 2*L[2]/dem
      mu= L[1]-sigma*( digamma((-eta-1)/h)
                      -digamma(1) ) -sigma*log(-h)
    }else{  
     num= 2*L[2]*k* (-h)^(k+1)
     dem= (eta+2)* ( (eta+1)*beta(k+1,(-eta-1)/h -k)
                    -(eta+2)*beta(k+1,(-eta-2)/h -k) )

     sigma = num/dem
     fac= 1- ( (eta+1)*beta(k+1,(-eta-1)/h -k) / ((-h)^(k+1)) )
     mu= L[1]-fac*sigma/k
    }

  }else if(h >0){
    
    if(abs(k) <= small){
      dem=(eta+2)*( digamma((eta+2)/h +1)
                   -digamma((eta+1)/h +1) ) 
      sigma= 2*L[2]/dem
      mu= L[1]-sigma*( digamma((eta+1)/h +1)
                      -digamma(1) ) - sigma*log(h)
    }else{
      num=2*L[2]*k* (h)^(k+1)
      dem= (eta+2)*( (eta+1)*beta(k+1,(eta+1)/h)
                    -(eta+2)*beta(k+1,(eta+2)/h) )
      sigma = num/dem
      fac= 1- ( (eta+1)*beta(k+1,(eta+1)/h) / (h)^(k+1) )
      mu= L[1]-fac*sigma/k
    }
  } # end if
  
  para=c(mu, sigma, k, h)
  
  tau.diff= abs(lhmom.kap(para,eta)$ratios[4]-t4)
  
  t34.diff= ( abs(lhmom.kap(para,eta)$ratios[3]-t3)
              + tau.diff )

  return(list(model = "k3d.kfix", para = para, eta=eta, 
              kfix=kfix, type="kap",
              tau.diff=tau.diff, ifail=ifail,
              precision=fvec, t34.diff=t34.diff,
              source="lh.park3k"))
}

#----------------------------------------------------
obj.lhk3k <- function(x, t3=t3, eta=eta,
                      kfix=kfix, hlow=hlow) { 
  
  #  h=x
  if ( x < hlow | x >5) return(10^6 )   # k for hosking style
  if ( x < 0 & (x*kfix) <= -eta-1) return(10^6) 
  #  if( k+0.725*h <= -1) return(c(10^6,10^6) ) 
  
  tau=rep(NA,4)
  tau[3:4]= lhmom.kap(c(0,1,kfix,x),eta)$ratios[3:4]

  return(tau[3]-t3)  
}
#---------------------------------------------------- 