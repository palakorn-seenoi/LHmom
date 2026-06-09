#------------------------------------------------------------------
# Estimating LH-me for Gumbel
# lh.pargum(data, eta=2)
#------------------------------------------------------------------
#' LH-moments parameters estimation for Gumbel distribution
#'
#' Function to calculate LH-moments parameters estimation for Gumbel distribution
#' @export
lh.pargum = function(data,eta=1){

  L= lhmoms(data, eta=eta)$lambdas[1:2]

  sigma= 2*L[2]/( (eta+2)*( log(eta+2)-log(eta+1) ) )
  mu=    L[1]-sigma*(0.5772 + log(eta+1))

  para=c(mu,sigma,0)
  names(para) <- c("xi","alpha","k")

  return( list(model="gum", para=para, type="gev",
               eta=eta, source="lh.pargum"
  ) )
}
