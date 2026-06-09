#---------------------------------------------
#' Initial parameters
#'
#' Function to calculate Initial parameters for LH-moments
#' @importFrom lmomco  lmoms
#' @export
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
