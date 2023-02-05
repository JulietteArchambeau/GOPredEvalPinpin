
# Function to identifiy outliers based on their loadings along significant RDA axes
# =================================================================================

# Forester et al. (2018)

# https://popgen.nescent.org/2018-03-27_RDA_GEA.html

detectoutliers <- function(x,z){                   # x = vector of loadings and z = number of standard deviations
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
