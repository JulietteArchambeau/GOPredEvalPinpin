# Evaluating genomic offset predictions in maritime pine


***

## REPORTS

The code (`.qmd` and `Rmd` files) used to generate the following reports can be found in the folder `/reports`.

-   [0_FormattingPopulationCoordinatesElevationClimateData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/0_FormattingPopulationCoordinatesElevationClimateData.html) Checking population information (coordinates and elevation data) from different sources - Extracting climatic data with ClimateDT - Calculating the average of the climatic variables across time periods of interest.

-   [1_FormattingGenomicData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/1_FormattingGenomicData.html) Formatting and filtering of the genomic data and imputation of missing data.

-   [2_CommonGardenData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/2_CommonGardenData.html). Extracting climatic data from ClimateDT and calculating the mean climate in each common garden between the planting date and the measurement date.

-   [3_CheckingPastFutureClimatesPopulationLocations.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/3_CheckingPastFutureClimatesPopulationLocations.html) Checking past and future climatic values at the population locations.

-   [4_ReduncancyAnalysis_VariancePartionning_IdentificationCandidateSNPs.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/4_ReduncancyAnalysis_VariancePartionning_IdentificationCandidateSNPs.html) Selection of the climatic variables, variance partitioning and identification of the candidate SNPs using Redundancy analysis (RDA) (approach developed in [Capblancq et al. 2018](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12906)).
    
    *  <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/PCAplot.pdf" target="_blank">PCAplot.pdf</a> Principal component analysis performed on the population allele frequencies. The first three axes of the PCA are used to account for the population structure in the RDA analysis.
 
    *  <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/RDAsummary.pdf" target="_blank">RDAsummary.pdf</a> Summary statistics of the RDA models.
    
    *  <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/RDAplots.pdf" target="_blank">RDAplots.pdf</a> RDA plots colored by gene pool.
    
    *  <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/RDAplots_outliers_1.pdf" target="_blank">RDAplots_outliers_1.pdf</a> RDA plots with outliers following [Forester et al. (2018)](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14584?casa_token=IOrVgFSER0gAAAAA%3AsOlFDnBLnWtTdC-R6vi5pZiRwuzpP4GQyr8H9hVpVqxW0_3RXOV6bznLQx9deVCrYv80LokfqFvaGeY) (and the [associated vignette](https://popgen.nescent.org/2018-03-27_RDA_GEA.html)).
    
    *  <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/RDAplots_outliers_2.pdf" target="_blank">RDAplots_outliers_2.pdf</a> RDA plots with outliers and  Manhattan plots following [Capblancq and Forester (2021)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13722) (and the [associated Github repository](https://github.com/Capblancq/RDA-landscape-genomics)).

-   [5_GradientForest_IdentificationCandidateSNPs.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/5_GradientForest_IdentificationCandidateSNPs.html) Identification of candidate SNPs with the Gradient Forest (GF) algorithm, using either raw allele frequencies (GF-raw) or allele frequencies after correction for population relatedness (GF-X), as described in [Fitzpatrick et al. 2021](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13374) and [Capblancq et al. 2023](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.18465).

-   [6_BaypassAnalysis_IdentificationCandidateSNPs.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/6_BaypassAnalysis_IdentificationCandidateSNPs.html) Identification of candidate SNPs with [BayPass](https://www1.montpellier.inra.fr/CBGP/software/baypass/index.html) (approach described in [Gautier 2015](https://academic.oup.com/genetics/article/201/4/1555/5930067?login=true)).

-   [7_LFMM_IdentificationCandidateSNPs_GenomicOffsetPredictions.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/7_LFMM_IdentificationCandidateSNPs_GenomicOffsetPredictions.html) Identification of candidate SNPs and estimation of the genetic gap (i.e. genomic offset) with the `LEA` R package ([Gain & Francois 2021](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13366)), which uses the latent factor mixed model (LFMM) approach ([Frichot et al. 2013](https://academic.oup.com/mbe/article/30/7/1687/972098?login=true); [Cayes et al. 2019](https://academic.oup.com/mbe/article/36/4/852/5290100?login=true)). Mathematical and theoretical details of the approach in [Gain et al. (2023)](https://www.biorxiv.org/content/10.1101/2023.01.02.522469v3).

-   [8_GeneratingSNPsets](https://juliettearchambeau.github.io/GOPredEvalPinpin/8_GeneratingSNPsets) Identifying the common candidates across the different gene-environment association (GEA) methods, looking at their genomic position and generating a set of control SNPs.

-   [9_GeneralizedDissimilarityModelling_GenomicOffsetPredictions.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/9_GeneralizedDissimilarityModelling_GenomicOffsetPredictions.html) Genomic offset predictions with the Generalized Dissimilarity Modelling (GDM) approach, as described in [Fitzpatrick & Keller 2015](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12376), [Mokany et al. 2022](https://onlinelibrary.wiley.com/doi/full/10.1111/geb.13459) and the [GDM website](https://mfitzpatrick.al.umces.edu/gdm/).

-   [10_GradientForest_GenomicOffsetPredictions.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/10_GradientForest_GenomicOffsetPredictions.html) Genomic offset predictions with the Gradient Forest (GF) algorithm, approach described in [Fitzpatrick & Keller 2015](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12376) and [Gougherty et al. 2021](https://www.nature.com/articles/s41558-020-00968-6).

    * <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/GFplots_cand.pdf" target="_blank">GFplots_cand.pdf</a> GF plots for the set of candidate SNPs identified by at least two GEA methods among RDA, pRDA, GF, LFMM and BayPass.
    
    * <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/GFplots_cand_corrected.pdf" target="_blank">GFplots_cand_corrected.pdf</a> GF plots for the set of candidate SNPs identified by at least two GEA methods correcting for population structure (i.e pRDA, LFMM, BayPass).
    
    * <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/GFplots_control.pdf" target="_blank">GFplots_control.pdf</a> GF plots for the set of control SNPs.
    
-   [11_RedundancyAnalysis_GenomicOffsetPredictions.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/11_RedundancyAnalysis_GenomicOffsetPredictions.html) Predicting the genomic offset with Redundancy Analysis (RDA).

-   [12_ComparingGenomicOffsetPredictions.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/12_ComparingGenomicOffsetPredictions.html) Comparing the genomic offset predictions across the different methods.

-   [13_ValidationNFI.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/13_ValidationNFI.html) Evaluating the genomic offset predictions with mortality rates in natural populations from the National Forest Inventory plots of France and Spain. The genomic offset predictions evaluated come from the four different methods (GF, GDM, RDA and LFMM) and are based on either the candidate or the control SNPs (and also all SNPs for LFMM).

-   [14_ValidationNFI_ModelComparison.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/14_ValidationNFI_ModelComparison.html) Building and evaluating the accuracy of the Bayesian models used to estimate the relationship between genomic offset predictions and mortality rates in the National Forest Inventory plots.

-   [15_ValidationCommonGardens.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/15_ValidationCommonGardens.html) Evaluating the genomic offset predictions with mortality and height data from five clonal common gardens (CLONAPIN network).
