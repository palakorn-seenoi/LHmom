library(lmomco)
library(nleqslv)

#---------- LH.par.k4d: theoretical LH moments ------
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
  
  xi=para[1]
  alpha = para[2]
  k=para[3]
  h=para[4]
  
  small=1e-5
  if(abs(k) < small) k=sign(k)*small
  if(abs(h) < small) h=sign(h)*small
  
  lambdas= ratios= rep(NA, nmom)
  B=rep(NA,8)
  
  if(k <= -1) stop("k should be > -1 in lhmom.k4d")
  
  if(h < 0){
    if( (k*h) <= -eta-1){
      stop("k*h should be > -eta-1 in lhmom.k4d")
      
    }else{
      for(ieta in 0:7){
        
        num= (ieta+1) #*gamma(k+1)*gamma((-ieta-1)/h -k)
        dem= ((-h)^(k+1))#*gamma((-ieta-1)/h +1)
        
        bb= num*beta(k+1,(-ieta-1)/h -k)/dem
        B[ieta+1]= xi + (alpha/k)*(1-bb)
      }
    } #end if k
    
  }else if(h >0){
    for(ieta in 0:7){
      num= (ieta+1)#*gamma(k+1)*gamma((ieta+1)/h )
      dem= (h^(k+1))#*gamma((ieta+1)/h +k+1)
      
      bb= num*beta(k+1,(ieta+1)/h)/dem
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
  
  z$lambdas = lam;  z$ratios = tau;  z$eta = eta
  z$type="kap"
  z
}
#---------------------------------------------------------
lh.parkap = function(data, eta=1, snap.tau4=FALSE,
                     nudge.tau4=1e-6, hlow= NULL,
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
        retrun(z)}
    }

    z$type="kap";  z$eta=0
    z$t34.diff= sum(abs(lhmom.kap(z$para,eta)$ratios[3:4]
                    -c(t3,t4)) )
    z$source="lh.park4d"
    z$snap.tau4 = snap.tau4; z$nudge.tau4=nudge.tau4

    z$ifailtext=NULL
    if(z$ifail== -1){
      z$ifailtext="para estim under snap.tau4 & nudge.tau4"
    }else if(z$ifail==0){
      z$ifailtext="Successful parameter estimation."
    }
    return(z)
  } # end if eta=0
  
  init=initkh.k4d(data, ntry=ntry)
  
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
    
#    cat("itry, termcd=", itry, mysol$termcd,"\n")
    
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
       z2= lh.park3k(data, eta,kfix=sol[1])
       
       if(z2$ifail==0) {
         ifail=1
         para=z2$para
         break
       }else{ 
         itry = itry + 1
         ftol= itry*ftolz
         xtol= itry*xtolz

         if (itry > ntry) { ifail =5
         break }
       }
        
      }else if(abs(fvec[1]) > ftol*10 & abs(fvec[2]) <= ftol*10){
       z2=list()
       z2= lh.park3h(data, eta,hfix=sol[2])
       
       if(z2$ifail==0) {
         ifail=1
         para=z2$para
         break
       }else{
         itry = itry + 1
         ftol= itry*ftolz
         xtol= itry*xtolz

         if (itry > ntry) { ifail =5
         break }
       }
       
      }else if(abs(fvec[1]) <= ftol*10 & abs(fvec[2]) <= ftol*10){
        ifail=0
        para[3:4] = mysol$x
        break
        
      }else{

        itry = itry + 1
        ftol= itry*ftolz
        xtol= itry*xtolz

        if (itry > ntry) { ifail =5
        break }
      }

    }else if(mysol$termcd >= 3){
      itry = itry + 1
      ftol= itry*ftolz
      xtol= itry*xtolz

      if (itry > ntry) { ifail =5
        break }
    } #end if sum
  } #end while
  }) # end trycatch
  
  if (ifail ==5) {
      # cat("failure to solve LHme for K4D, eta=",eta,"\n")
      # zg=lh.pargev(data,eta)
      
      z=list()
      # z$para=c(zg$para, -small)
      
#      z= lh.park3h(data,eta, hfix= hlow+subt)
      z$para=rep(NA,4)
      z$ifail =ifail; z$type="kap"
      # z$precision =mysol$fvec
      z$eta=eta
      # z$t34.diff= sum(abs(lhmom.kap(z$para,eta)$ratios[3:4]
      #               -c(t3,t4)) )
      z$snap.tau4 = snap.tau4; z$nudge.tau4=subt
      z$ifailtext="failue for k4d"
      return(z)
  }

  L= samlh$lambdas[1:2]
  k= para[3]
  h= para[4]
  
#  if(h < -1) h=-0.99999 -------------------------
  if(h < hlow) h= hlow +small
  if(k < -1) k=-0.99999

if(ifail==0){
 if(para[4] <0){
  num=2*L[2]*k* (-h)^(k+1)
  dem= (eta+2)* ( (eta+1)*beta(k+1,(-eta-1)/h -k)
                 -(eta+2)*beta(k+1,(-eta-2)/h -k) )
                # gamma(k+1)*( (eta+1)*G[1]-(eta+2)*G[2] ) 
  sigma = num/dem
  fac= 1- ( (eta+1)*beta(k+1,(-eta-1)/h -k) / ((-h)^(k+1)) )
  mu= L[1]-fac*sigma/k
  
 }else if(para[4] >0){
  num=2*L[2]*k* (h)^(k+1)
  dem= (eta+2)*( (eta+1)*beta(k+1,(eta+1)/h)
                -(eta+2)*beta(k+1,(eta+2)/h) )
  sigma = num/dem
  fac= 1- ( (eta+1)*beta(k+1,(eta+1)/h) / (h)^(k+1) )
  mu= L[1]-fac*sigma/k
 } # end if
 para= c(mu, sigma, k, h)
}

tau.diff= sum(abs(lhmom.kap(para,eta)$ratios[3:4]
              -c(t3,t4)))

# if(ifail==6){
#   ifailtext="para estim under snap.tau & nudge.tau4"
# }else 
if(ifail==5){
  ifailtext="failue for k4d"
}else if(ifail==0){
  ifailtext="Successful parameter estimation."
}else if(ifail==1){
  ifailtext=paste("x within xtol, but fcn value is > ftol.",
          "para are obtained using lh.park3k or lh.park3h")
}

return(list(type = "kap", para = para, eta=eta, 
            t34.diff=tau.diff, ifail=ifail,
            precision=fvec, snap.tau4=snap.tau4,
            nudge.tau4=subt, source="lh.parkap",
            ifailtext=ifailtext))
}

# x=as.vector(init[2, 1:2])
#----------------------------------------------------
obj.lhk4d <- function(x, t3=t3,t4=t4, eta=eta,
                      hlow=hlow) { 
  
  k=x[1]
  h=x[2]

  if ( k <= -1 ) return(rep(10^6,2) )   # k for hosking style
  if (h < 0 & (k*h) <= -eta-1) return(c(10^6,10^6) ) 

  if(h < hlow | h > 5 ) return(c(10^6,10^6) )
  
#  if( k+0.725*h <= -1) return(c(10^6,10^6) ) 

  tau=rep(NA,4)
  tau[3:4]= lhmom.kap(c(0,1,k,h),eta)$ratios[3:4]
  
  if( any(is.na(tau[3:4])) ) return(rep(10^6,2))
  
  # if( tau[4] >= (5*(tau[3]^2)+1)/6) { 
  #   if(subt==0) { return(c(10^6,10^6) )
  #   }else if(subt != 0){
  #     tau[4]= (5*(tau[3]^2)+1)/6 - subt  }
  # }

  return(c(tau[3]-t3,tau[4]-t4))
}
#---------------------------------------------------- 
#--------------------------------------------------------  
initkh.k4d = function(data, ntry=5){
  
  init= matrix(NA,max(5,ntry),2)
  lmom0=lmoms(data)

  init[1,]= parkap(lmom0, snap.tau4=TRUE,
                   nudge.tau4=1e-3)$para[3:4]
  if(any(is.na(init[1,1:2]))) {
    init[1,]= c(pargev(lmom0)$para[3],-0.05) }
    
  init[2,]= c(-0.25,-0.75)
  init[3,]= park3h(lmoms(data),hfix=-2)$para[3:4]
  init[4,]= park3h(lmoms(data),hfix=1.5)$para[3:4]
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
# #------------------------------------------------------
# lhmom.gpa = function(para=NULL,eta=1){
#   lhmom.kap(c(para,1),eta)
# }
# lhmom.glo = function(para=NULL,eta=1){
#   lhmom.kap(c(para,-1),eta)
# }
