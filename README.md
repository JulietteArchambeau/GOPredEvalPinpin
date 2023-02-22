# Evaluating genomic offset predictions in maritime pine


***

## REPORTS

The code (`.qmd` and `Rmd` files) used to generate the following reports can be found in the folder `/reports`.

-   [0_FormattingPopulationCoordinatesElevationClimateData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/0_FormattingPopulationCoordinatesElevationClimateData.html) Checking population information (coordinates and elevation data) from different sources - Extracting climatic data with ClimateDT.

-   [1_FormattingGenomicData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/1_FormattingGenomicData.html) Formatting and filtering of the genomic data and imputation of missing data.

-   [2_ReduncancyAnalysis_part1.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/2_ReduncancyAnalysis_part1.html) Selection of the climatic variables, variance partitioning and identification of the candidate SNPs using RDA (Redundancy analysis)
    
    Associated figures and documents:
 
    *  <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/RDAsummary.pdf" target="_blank">RDAsummary.pdf</a> Summary statistics of the RDA models.
    
    *  <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/RDAplots.pdf" target="_blank">RDAplots.pdf</a> RDA plots colored by gene pool.
    
    *  <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/RDAplots_outliers_1.pdf" target="_blank">RDAplots_outliers_1.pdf</a> RDA plots with outliers following [Forester et al. (2018)](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14584?casa_token=IOrVgFSER0gAAAAA%3AsOlFDnBLnWtTdC-R6vi5pZiRwuzpP4GQyr8H9hVpVqxW0_3RXOV6bznLQx9deVCrYv80LokfqFvaGeY) (and the [associated vignette](https://popgen.nescent.org/2018-03-27_RDA_GEA.html)).
    
    *  <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/RDAplots_outliers_2.pdf" target="_blank">RDAplots_outliers_2.pdf</a> RDA plots with outliers and  Manhattan plots following [Capblancq and Forester (2021)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13722) (and the [associated Github repository](https://github.com/Capblancq/RDA-landscape-genomics)).

-   [3_GradientForest_part1.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/3_GradientForest_part1.html) Identification of candidate SNPs with the Gradient Forest algorithm.


