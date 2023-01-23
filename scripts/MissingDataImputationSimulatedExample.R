# ===========================================
# Comparing two ways of imputing missing data
# ===========================================

# 23/01/2023


# We want to impute missing data on the basis on the main gene pool of each individual.


library(tidyverse)
library(gtools)
library(tictoc)

set.seed(485)


# Data simulation
# ===============

# Not realistic at all, but will do the job :D

gp.nb <- 6 # Number of gene pools
ind.nb.per.gp <- 10 # Number of individuals in each gene pool
nb.snp <- 100 # number of genetic markers (ie SNPs)


# Function to generate a vector of allele counts with some NAs
GenerateAlleleCounts <- function(x){
  
  vec <- sapply(1:gp.nb, function(x) c(0,1,2) %>% 
                  sample(ind.nb.per.gp, replace = TRUE, prob=as.vector(rdirichlet(1,c(1/3,1/3,1/3))))) %>% 
    as.vector() 
  
  nb.na <- rbinom(n=1,size=60,prob=0.05)
  id.na <- sample(1:(gp.nb*ind.nb.per.gp),nb.na,replace=F)
  vec[id.na] <- NA
  
  return(vec)
  
}

df <- sapply(1:nb.snp, GenerateAlleleCounts) %>% 
  as.data.frame() %>% 
  set_colnames(paste0("snp",1:nb.snp))


df <- data.frame(individual=paste0("i",1:(gp.nb*10)), 
                 gene.pool=sapply(1:gp.nb, function(x) rep(paste0("gp",x),ind.nb.per.gp)) %>% as.vector()) %>% 
  bind_cols(df)



# Using a loop for to impute missing data
# =======================================

dfloop <- df
  
tic("Loop")
for(i in unique(dfloop$gene.pool)){
    subset <- dfloop[dfloop$gene.pool==i,]
    subset <- apply(subset[,3:ncol(subset)], 2, function(x) replace(x, is.na(x), as.numeric(names(which.max(table(x))))))
    dfloop[dfloop$gene.pool==i,3:ncol(dfloop)] <- subset
}
toc()



# Using tidyverse functions for to impute missing data
# ====================================================

tic("Tidyverse")
dftidy <- df %>% 
  group_by(gene.pool) %>% 
  group_split() %>% 
  map_dfr(function(gene.pool){
    gene.pool %>% 
      mutate_if(is.numeric, ~replace_na(.,as.numeric(names(which.max(table(.))))))
  })
toc()



# We check that the two datasets produced by the two functions are identical.
identical(as.data.frame(dftidy),dfloop) # that's ok :)


