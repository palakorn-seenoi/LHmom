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

# lhmom.gno(para=c(100,10,-0.3), eta=0)
# ==============================================================
# Function to calculate Theoretical LH-moments 
# ==============================================================

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
  lambdas[2] <- alpha*(eta+2)/(2*k) * (B1 -B2) 

  lambdas[3] <- alpha*(eta+3)/(6*k) * ( -(eta+2)*B1 + 2*(eta+3)*B2
                                        - (eta+4)*B3 )
  
  # lambda4_eta
  lambdas[4] <- alpha*(eta+4)/(24*k) * ( (eta+2)*(eta+3)*B1
                                          - 3*(eta+3)*(eta+4)*B2
                                          + 3*(eta+4)*(eta+5)*B3
                                          -   (eta+5)*(eta+6)*B4 )
  
  ratios[2]=lambdas[2]/lambdas[1]
  for (j in 3:nmom) ratios[j]=lambdas[j]/lambdas[2]
  
  z=list()
  z$lambdas = lambdas;  z$ratios = ratios;  z$eta = eta
  z$type="gno"
  z
}

# lh.pargno(data,eta=4)
# =============================================================
# LHme for gno
# =============================================================
lh.pargno <- function(data, eta=1, opt=FALSE, ntry=5){ 

  z=list()
  if(length(eta) != 1) stop("eta should be a scalar")
  if( eta %% 1 != 0 | eta < 0 | eta >= 5) {
    stop("eta must be a nonnegative integer less than 5")}
  
  samlh=lhmoms(data, eta=eta)
  t3= samlh$ratios[3]
  t4= samlh$ratios[4]
  
  if(eta==0){
    z= pargno(lmoms(data))
    z$eta=0
    z$tau.diff= abs(lhmom.gno(z$para,eta)$ratios[4]-t4)
    z$t34.diff= ( abs(lhmom.gno(z$para,eta)$ratios[3]-t3)
                + z$tau.diff )
    return(z)
  }
  
if(opt != TRUE){     # obtain LH-me without numerical solution
 
    aa=matrix(NA,5,4)
    # aa[1,]= c(0,	-2.0462,	0,	-0.5016)
    aa[2,]= c(0.2825,	-2.5672,	-0.1521,	-0.6662)
    aa[3,]=	c(0.4720,	-2.9385,	-0.2351,	-0.9060)
    aa[4,]=	c(0.6142,	-3.2326,	-0.2694,	-1.1471)
    aa[5,]=	c(0.7276,	-3.4776,	-0.2746,	-1.3719)
    
    jeta=eta+1
    k= aa[jeta,1]+ aa[jeta,2]*t3 + aa[jeta,3]*t3^2 +aa[jeta,4]*t3^3
    ifail=0
    
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
  tau.diff=abs(lhmom.gno(para,eta)$ratios[4]-t4)
  
  t34.diff= ( abs(lhmom.gno(para,eta)$ratios[3]-t3)
             + tau.diff )
  
  return(list(type = "gno", para = para, eta=eta, 
              tau.diff=tau.diff, source="lh.pargno",
              t34.diff=t34.diff, ifail=ifail))
}

#----------- sample running ---------------------------------------
# eta=0
# pargno(lmoms(data))
# 
# for (eta in 0:4){
#   lhold= lh.pargno(data,eta=eta, opt=FALSE)
#   lhnew= lh.pargno(data,eta=eta, opt=TRUE)
#   cat("eta, para.old,para.new=",eta,lhold$para,lhnew$para,"\n")
# }



#---------------------------------------------
# Seenoi

data_rainfall <- read.csv("D:/SoftwareX/rainfall_max_thailand/max_rainfall_khonkaen.csv")
data <- data_rainfall$rainfall
data
eta=1
pargno_est <- pargno(lmoms(data))$para

ks.test(data, ) 
 
for (eta in 0:4){
   lhold= lh.pargno(data,eta=eta, opt=FALSE)
   lhnew= lh.pargno(data,eta=eta, opt=TRUE)
   cat("eta, para.old,para.new=",eta,lhold$para,lhnew$para,"\n")
}
