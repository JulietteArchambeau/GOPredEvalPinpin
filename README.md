# Evaluating genomic offset predictions in maritime pine


***

## REPORTS

The code (`.qmd` and `Rmd` files) used to generate the following reports can be found in the folder `/reports`.

-   [0_FormattingPopulationCoordinatesElevationClimateData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/0_FormattingPopulationCoordinatesElevationClimateData.html) Checking population information (coordinates and elevation data) from different sources - Extracting climatic data with the [Climate Downscaling Tool (ClimateDT)](https://www.ibbr.cnr.it/climate-dt/) - Calculating the average of the climatic variables across time periods of interest.

-   [1_FormattingGenomicData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/1_FormattingGenomicData.html) Formatting and filtering of the genomic data, checking SNPs position on the genome, imputation of missing data.

-   [2_CommonGardenData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/2_CommonGardenData.html). Extracting climatic data from [ClimateDT](https://www.ibbr.cnr.it/climate-dt/) at the location of the common gardens and calculating the mean climate in each common garden between the planting date and the measurement date.

-   [3_CheckingPastFutureClimatesPopulationLocations.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/3_CheckingPastFutureClimatesPopulationLocations.html) Checking differences in past and future climatic values at the locations of the populations.

-   [4_ReduncancyAnalysis_VariancePartionning_IdentificationCandidateSNPs.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/4_ReduncancyAnalysis_VariancePartionning_IdentificationCandidateSNPs.html) Selection of the climatic variables, variance partitioning and identification of the candidate SNPs using Redundancy analysis (RDA) (approach developed in [Capblancq et al. 2018](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12906) and [Capblancq and Forester 2021](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13722)).
    
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

-   [12_ComparingGenomicOffsetPredictions.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/12_ComparingGenomicOffsetPredictions.html) Comparing the genomic offset predictions across the different methods, SNPs sets and Global Climate Models (GCMs).

-   [13_ValidationNFI.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/13_ValidationNFI.html) Evaluating the genomic offset predictions with mortality rates in natural populations from the National Forest Inventory plots of France and Spain. Genomic offset predictions come from the four different methods (GF, GDM, RDA and LFMM) and are based on either the candidate or the control SNPs (and also all SNPs for LFMM).

-   [14_ValidationNFI_ModelComparison.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/14_ValidationNFI_ModelComparison.html) Building and evaluating the accuracy of the Bayesian models used to estimate the relationship between genomic offset predictions and mortality rates in the National Forest Inventory plots.

-   [15_ValidationCommonGardens.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/15_ValidationCommonGardens.html) Evaluating the genomic offset predictions with mortality and height data from five clonal common gardens (CLONAPIN network).

-   [16_PCAplotIDpops.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/16_PCAplotIDpops.html) Generating figures for the Supplementary Information based on the PCA of the control and candidate SNPs: screeplots and PCA plots with the ALT and ARM populations highlighted.

-   [RPackageCitations.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/RPackageCitations.html) Citations of the R packages used in the present study.

***

## Description of the DRYAD repository

Datasets associated with the paper were deposited in the following DRYAD repository:


### Genomic data

#### Population genetic structure

**Dataset `PopulationStructureCorrea2015.csv`**

Proportion of assignement to the six gene pools identified in [Jaramillo-Correa et al. (2015)](https://academic.oup.com/genetics/article/199/3/793/5935834?login=false). This dataset contains 523 genotypes but only 454 genotypes were used in the present study because the other genotypes had no genomic data or too much missing data.


Meaning of the columns:

  1. `clon`: clone (i.e. genotype)
  2. `pop`: population (i.e. provenance)
  3. `Q1`: proportion of assignment to the northern African (NA) gene pool for each clone.
  4. `Q2`: proportion of assignment to the Corsican (C) gene pool for each clone.
  5. `Q3`: proportion of assignment to the central Spain (CS) gene pool for each clone.
  6. `Q4`: proportion of assignment to the French Atlantic (FA) gene pool for each clone.
  7. `Q5`: proportion of assignment to the Iberian Atlantic (IA) gene pool for each clone.
  8. `Q6`: proportion of assignment to the south-eastern Spain (SES) gene pool for each clone.
  9. `max.Q`: ID of the main gene pool of each clone.
  10. `main_gp`: full name of the main gene pool of each clone.
  11. `color_main_gp`: color attributed to the main gene pool of each clone in the figures.


&nbsp;

#### Raw genomic data

**Dataset `RawGenomicData.csv`**

Genomic data before data filtering and formatting. NAs are indicated with `---`.

Meaning of the columns:

  1. `clone`: clone ID.
  2. `assay`: assay in which the clone was genotyped, either the Infinium assay (`only_Inf`), the Axiom assay (`only_Affx`) or both assays (`both_Inf_Affx`).
  3. `snp_1` to `snp_14016`:  genotype for each of the 14,016 SNPs.


&nbsp;


#### Filtered and formatted genomic data

For data filtering and formatting, see report [1_FormattingGenomicData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/1_FormattingGenomicData.html).

**Dataset `FormattedFilteredGenomicData_AlleleCounts_withmaf.csv`**

Genomic information of each genotype (coded as 0, 1 or 2) including minor allele frequencies (MAF), i.e. MAF < 1%. SNPs in rows and genotypes in columns.


**Dataset `FormattedFilteredGenomicData_AlleleCounts_withoutmaf.csv`**

Genomic information of each genotype (coded as 0, 1 or 2) without MAF. SNPs in rows and genotypes in columns.


**Dataset `FormattedFilteredGenomicData_AlleleFrequencies_withoutmaf.csv`**

Allele frequencies of the populations without MAF. SNPs in columns and populations in rows.

&nbsp;

#### Imputed genomic data

For genomic data imputation, see section 6 of the report [1_FormattingGenomicData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/1_FormattingGenomicData.html).

**Dataset `ImputedGenomicData_AlleleCounts_withmaf.csv`**

Imputed genomic data of each genotype (coded as 0, 1 or 2) including MAF. SNPs in rows and genotypes in columns.

**Dataset `ImputedGenomicData_AlleleCounts_withmaf.csv`**

Imputed genomic data of each genotype (coded as 0, 1 or 2) without MAF. SNPs in rows and genotypes in columns.

**Dataset `ImputedGenomicData_AlleleFrequencies_withmaf.csv`**

Imputed allele frequencies of the populations with MAF. SNPs in columns and populations in rows.

**Dataset `ImputedGenomicData_AlleleFrequencies_withoutmaf.csv`**

Imputed allele frequencies of the populations without MAF. SNPs in columns and populations in rows.

&nbsp;

#### Genomic position and other SNP IDs


**Dataset `ListSNPs_withmaf.csv`**

Supplementary information about all the SNPs used in the present paper (including SNPs with MAF < 1%). SNPs are in rows. Meaning of the columns:

  1. `snp_ID`: SNP ID used in the present study.
  2. `original_ID`: Original ID of the SNP.
  3. `affx_ID`: SNP ID in the Axiom assay.
  4. `infinium_ID`: SNP ID in the Infinium assay.
  5. `scaffold/contig`: contig on which the SNP is located. The term scaffold is also used as some SNPs were obtained from the alignment of NGS short-reads from a pseudoreference genome in *Pinus pinaster* that were called scaffolds (even though there are not really scaffolds).
  6. `genome_position`: position of the SNP on the scaffold/contig (in bp).
  7. `annotation`: SNP name (including possible alternative names).

**Dataset `ListSNPs_withoutmaf.csv`**

Supplementary information about the SNPs with MAF > 1%. SNPs are in rows. See dataset `ListSNPs_withmaf.csv` for the meaning of the columns.

&nbsp;

### Climatic data

Climatic data was extracted with the [Climate Downscaling Tool (ClimateDT)](https://www.ibbr.cnr.it/climate-dt/). For the full name and metadata of the climatic variables, see https://www.ibbr.cnr.it/climate-dt/?action=fldlist.


**Dataset `ClimateDT_CommonGardens.csv`**

Climatic data (from 2010 to 2018) at the location of the five common gardens: Asturias (Spain), Bordeaux (France), Cáceres (Spain), Madrid (Spain) and Fundão (Portugal). See report [2_CommonGardenData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/2_CommonGardenData.html).


**Dataset `ClimateDT_Populations_PastClimates.csv`**

Climatic data (from 1901 to 2021) at the location of the 34 populations. See report [0_FormattingPopulationCoordinatesElevationClimateData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/0_FormattingPopulationCoordinatesElevationClimateData.html).


**Dataset `ClimateDT_Populations_FutureClimates.csv`**

Predictions of future climatic data  at the location of the 34 populations for the period 2041-2070, under the scenario SSP3.7-0 and for five Global Climate Models (GCMs): GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0 and UKESM1-0-LL (column `gcm`). See report [0_FormattingPopulationCoordinatesElevationClimateData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/0_FormattingPopulationCoordinatesElevationClimateData.html).


**Dataset `ClimateDT_NFIPlots_PastClimates.csv`**

Averaged climatic data  at the location of the NFI plots for the reference period 1901-1950. See report [13_ValidationNFI.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/13_ValidationNFI.html).


**Dataset `ClimateDT_NFIPlots_SurveyClimates.csv`**

Averaged climatic data  at the location of the NFI plots for the survey periods specific to each plot. See report [13_ValidationNFI.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/13_ValidationNFI.html).

&nbsp;

### National Forest Inventory (NFI) data

**Dataset `NFIdata_cleaned.csv`**

This dataset contains mortality data from the natural populations of the Spanish and French NFI plots in which maritime pines were recored. See report [13_ValidationNFI.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/13_ValidationNFI.html). Meaning of the columns:

  1. `plotcode`: code of the NFI plot.
  2. `longitude`: longitude of the NFI plot.
  3. `latitude`: latitude of the NFI plot.
  4. `country`: country of the NFI plot (ES = Spain, FR = France).
  5. `nb_years`: number of years between surveys. It corresponds to the number of years between the first and second inventory in Spain and is equal to 5 in the French inventory as mortality was estimated in the five years before the survey date.
  6. `first_survey`: year of the first survey in the Spanish NFI (`NA` for the French NFI)
  7. `second_survey`: year of the second survey in the Spanish NFI and the unique survey in the French NFI
  8. `nb_dead`: number of dead maritime pines in the NFI plot.
  9. `nb_tot`: total number of maritime pines in the NFI plot (including the dead ones).
  10. `basal_area`: basal area of all tree species in the NFI plot (proxy of the competition among trees).
  11. `prop_dead`: proportion of dead maritime pines in the NFI plot.
  12. `annual_prop_dead`: annual proportion of maritime pine mortality in the NFI plot.

&nbsp;

### Common garden data

See report [15_ValidationCommonGardens.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/15_ValidationCommonGardens.html)

**Dataset `CommonGardendata_cleaned.csv`**

This dataset contains mortality and height data from five clonal common gardens. Meaning of the columns:

  1. `cg`: name of the common garden.
  2. `block`: block in which the tree was planted.
  3. `pop`: population of the tree.
  4. `clon`: genotype of the tree.
  5. `tree`: tree ID.
  6. `AST_htmar14`: height measurements in March 2014 in Asturias (Spain) when the trees were 37 month-old (trees planted in February 2011).
  7. `BDX_htnov2018`: height measurements in November 2018 in Bordeaux (France) when the trees were 85 month-old (trees planted in October 2011).
  8. `CAC_htdec11`: height measurements in December 2011 in Cáceres (Spain) when the trees were 8 month-old (trees planted in April 2011).
  9. `MAD_htdec11`: height measurements in December 2011 in Madrid (Spain) when the trees were 13 month-old (trees planted in November 2010).
  10. `POR_htmay13`: height measurements in May 2013 in Fundão (Portugal) when the trees were 27 month-old (trees planted in February 2011).
  11. `AST_survmar14`: survival data (`0` for dead trees and `1` for survivors) in March 2014 in Asturias.
  12. `BDX_surv18`: survival data in November 2018 in Bordeaux. 
  13. `CAC_survdec11`: survival data in December 2011 in Cáceres.
  14. `MAD_survdec11`: survival data in December 2011 in Madrid.
  15. `POR_survmay13`: survival data in May 2013 in Fundão.
  

**Dataset `HeightIntercepts_Archambeauetal2023.csv`**

This dataset contains the population height intercepts calculated across the five common gardens in the model 1 of [Archambeau et al. (2022)](https://www.journals.uchicago.edu/doi/10.1086/720619). Meaning of the columns:

  1. `pop`: population ID.
  2. `height`: mean of the posterior distributions of the population varying intercepts.
  3. `std.error`: standard error of the mean of the posterior distributions of the population varying intercepts.
  4. and 5. `conf.low` and `conf.high`: credible intervals of the posterior distributions of the population varying intercepts.
