
#------------------------------------------------------------------
# Calculate sample LH-moments from data, using Wang's paper(1997)
# lhmoms(data,eta=seq(0,4)) or lhmoms(data,eta=2)
#------------------------------------------------------------------

#' Sample LH-moments from data
#'
#' Calculate sample LH-moments from sample data
#' @export
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
    z$eta=as.vector(eta)
    z$lambdas =as.matrix(t(lambda[,jeta]), nrow=1);
    z$ratios = as.matrix(t(ratio[,jeta]), nrow=1)
    rownames(z$lambdas) <- paste0("eta=",eta)
    rownames(z$ratios)  <- paste0("eta=",eta)

  }else{
    z$eta=eta
    z$lambdas = t(lambda)
    z$ratios = t(ratio)
    rownames(z$lambdas) <- paste0("eta=",1:leta)
    rownames(z$ratios)  <- paste0("eta=",1:leta)
  }
  colnames(z$lambdas) <- paste0("lhmom-",1:nmom)
  colnames(z$ratios)  <- paste0("lht-",1:nmom," ")
  z
}
