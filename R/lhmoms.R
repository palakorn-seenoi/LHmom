
#------------------------------------------------------------------
# Calculate sample LH-moments from data, using Wang's paper(1997)
# lhmoms(data,eta=seq(0,4)) or lhmoms(data,eta=2)
#------------------------------------------------------------------

#' Calculate Sample LH-moments
#'
#' This function computes the sample LH-moments and LH-moment ratios from a numeric
#' dataset. The estimation is based on the methodology proposed by Wang (1997).
#' The order parameter \code{eta} determines the weight assigned to larger observations.
#' When \code{eta = 0}, the calculation reduces to the ordinary sample L-moments.
#'
#' @param data A numeric vector of data values.
#' @param eta A non-negative integer (or a sequence of integers) between 0 and 4
#'   representing the order of the LH-moments.
#' @param nmom An integer specifying the maximum number of moments to compute (default is 5).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{eta}: The order(s) of the LH-moments.
#'   \item \code{lambdas}: A matrix of the calculated sample LH-moments.
#'   \item \code{ratios}: A matrix of the calculated sample LH-moment ratios.
#' }
#'
#' @references Wang, Q. J. (1997). Using higher order L-moments for regional
#' flood frequency analysis. \emph{Water Resources Research}, 33(12), 2841-2848.
#'
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

