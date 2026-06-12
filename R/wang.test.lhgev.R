
#----------------------------------------------------------------
# Perform Wang's GOF test for GEVD, using LH moments at WRR(1998)
# wang.test.lhgev(data)
#----------------------------------------------------------------

#' Wang's Goodness-of-Fit Test for the Generalized Extreme Value (GEV) Distribution
#'
#' This function performs a Goodness-of-Fit (GOF) test for the Generalized Extreme
#' Value (GEV) distribution using LH-moments, based on the methodology proposed
#' by Wang (1998). It calculates a Z-test statistic by comparing the sample
#' LH-kurtosis with the theoretical LH-kurtosis. The test is evaluated across
#' LH-moment orders (\code{eta}) from 0 to 4.
#'
#' @param data A numeric vector of data values.
#'
#' @return A data frame containing the GOF test results for each \code{eta} value
#'   (from 0 to 4), with the following columns:
#' \itemize{
#'   \item \code{eta}: The order of the LH-moments.
#'   \item \code{z.test}: The calculated Z-statistic for the GOF test.
#'   \item \code{cond.sigma}: The conditional standard deviation of the sample LH-kurtosis.
#'   \item \code{p.value}: The two-sided p-value corresponding to the Z-test statistic.
#' }
#'
#' @references Wang, Q. J. (1998). Approximate goodness-of-fit tests of fitted
#' generalized extreme value distributions using LH moments. \emph{Water Resources
#' Research}, 34(12), 3497-3502.
#'
#' @export
wang.test.lhgev =function(data){

  # tau = population L-moments
  # t = sample L-moment
  # n= sample size

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

  n=length(data)
  result=matrix(NA,5,4)

  for(ieta in 0:4){

    LH=  lh.pargev(data,eta=ieta)       # computing LHme for gevd
    tau= lhmom.gev(LH$para,eta=ieta)$ratios

    t= lhmoms(data,eta=ieta)$ratios    # sample L-moments and L-ratios

    # lh.pargev: computing LHme for gevd
    # lhmom.gev: population L-moments and L-ratios for gevd
    # lhmoms: sample L-moments and L-ratios for given eta

    jeta= ieta+1
    B= bb[jeta,1] +bb[jeta,2]*tau[3] +bb[jeta,3]*tau[3]^2
    B= B +bb[jeta,4]*tau[3]^3 +bb[jeta,5]*tau[3]^4

    C= cc[jeta,1] +cc[jeta,2]*tau[3] +cc[jeta,3]*tau[3]^2
    C= C +cc[jeta,4]*tau[3]^3 +cc[jeta,5]*tau[3]^4

    cond.sigma_t4 = sqrt(B/n + C/(n^2))

    result[jeta,1]= jeta-1
    result[jeta,2]= (t[4]-tau[4])/cond.sigma_t4
    result[jeta,3]= cond.sigma_t4
    result[jeta,4]= 2*pnorm( abs( result[jeta,2]),lower.tail=FALSE)
  }

  colnames(result) = c("eta","z.test", "cond.sigma","p.value")
  as.data.frame(result)
}
