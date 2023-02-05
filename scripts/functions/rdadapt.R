
# Function to conduct a RDA based genome scan
# ===========================================

# From Capblancq & Forester (2021)

# https://github.com/Capblancq/RDA-landscape-genomics/blob/main/src/rdadapt.R 

# required library: library(robust)

rdadapt <- function(rda,K)
{
  zscores<-rda$CCA$v[,1:as.numeric(K)]
  resscale <- apply(zscores, 2, scale)
  resmaha <- covRob(resscale, distance = TRUE, na.action= na.omit, estim="pairwiseGK")$dist
  lambda <- median(resmaha)/qchisq(0.5,df=K)
  reschi2test <- pchisq(resmaha/lambda,K,lower.tail=FALSE)
  qval <- qvalue(reschi2test)
  return(data.frame(pvalues=reschi2test, qvalues=qval$qvalues))
}
