
# Function to identifiy outliers based on their loadings along significant RDA axes
# =================================================================================

# from Forester et al. (2018)

# Forester BR, Lasky JR, Wagner HH, Urban DL (2018) Comparing methods for detecting multilocus adaptation with multivariate genotype-environment associations. Molecular Ecology 27(9):2215-2233.

# https://popgen.nescent.org/2018-03-27_RDA_GEA.html

detectoutliers <- function(x,z){             # x = vector of loadings and z = number of standard deviations
  lims <- mean(x) + c(-1, 1) * z * sd(x)     # find loadings +/-z sd from mean loading     
  x[x < lims[1] | x > lims[2]]               # locus names in these tails
}
