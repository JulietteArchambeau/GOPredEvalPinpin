# Data and code for the paper: 'Evaluating genomic offset predictions in a forest tree with high population genetic structure'

***

Juliette Archambeau<sup>1,2</sup>, Marta Benito-Garzón<sup>1</sup>, Marina de-Miguel<sup>1,3</sup>, Alexandre Changenet<sup>1</sup>, Francesca Bagnoli<sup>4</sup>, Frédéric Barraquand<sup>5</sup>, Maurizio Marchi<sup>4</sup>, Giovanni G. Vendramin<sup>4</sup>, Stephen Cavers<sup>2</sup>, Annika Perry<sup>2</sup> and Santiago C. González-Martínez<sup>1</sup>

**1** INRAE, Univ. Bordeaux, BIOGECO, F-33610 Cestas, France

**2** UK Centre for Ecology \& Hydrology, Bush Estate, Penicuik, United Kingdom

**3** EGFV, Univ. Bordeaux, Bordeaux Sciences Agro, INRAE, ISVV, F-33882, Villenave d'Ornon, France

**4** Institute of Biosciences and BioResources, National Research Council, 50019 Sesto Fiorentino, Italy

**5** CNRS, Institute of Mathematics of Bordeaux, F-33400 Talence, France

**Corresponding author:** Juliette Archambeau, juli.archambeau@gmail.com

***

**Paper abstract**

Predicting how tree populations will respond to climate change is an urgent societal concern. An increasingly popular way to make such predictions is the genomic offset (GO) approach, which uses current gene-environment associations to identify populations that may experience climate maladaptation in the near future. However, GO has strong limitations and, despite promising validation of its predictions using height data from common gardens, it still lacks broad empirical testing. Using maritime pine, a tree species from southwestern Europe and North Africa with a marked population genetic structure, we evaluated GO predictions from four methods, namely Gradient Forest (GF), Redundancy Analysis (RDA), latent factor mixed models (LFMM) and Generalised Dissimilarity Modeling (GDM). GO was predicted using 9,817 SNPs genotyped on 454 trees from 34 populations and was then validated with mortality data from National Forest Inventories and mortality and height data from five common gardens. We found high variability in GO predictions and validation. GO predictions with GDM and GF (and to a lesser extent RDA) based on candidate SNPs potentially involved in climate adaptation showed the strongest and most consistent associations with mortality rates in common gardens and NFI plots. We found almost no association between GO predictions and tree height in common gardens, most likely due to the overwhelming effect of population genetic structure on tree height. Our study demonstrates the imperative to validate GO predictions with a range of independent data sources before they can be used as informative and reliable metrics in conservation or management strategies.

***

## Data in the DRYAD repository

### Genomic data

#### Population genetic structure

**Dataset `PopulationStructureCorrea2015.csv`**

Proportion of assignement to the six gene pools identified in [Jaramillo-Correa et al. (2015)](https://academic.oup.com/genetics/article/199/3/793/5935834?login=false). This dataset contains 523 genotypes but only 454 genotypes were used in the present study because the other genotypes had no genomic data or too much missing data.


Meaning of the columns:

  1. `clon`: clone (i.e., genotype)
  2. `pop`: population (i.e., provenance)
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

For data filtering and formatting, see `1_FormattingGenomicData.qmd`.

**Datasets `FormattedFilteredGenomicData_AlleleCounts_withmaf.csv` and `FormattedFilteredGenomicData_AlleleCounts_withoutmaf.csv`**

Allele counts of each genotype (coded as 0, 1 or 2). SNPs in rows and genotypes in columns.

  1. `snp_ID`: SNP ID.
  2. to 455. `ALT10` to `VER9`:  Clone (i.e., genotype) ID.
  
`FormattedFilteredGenomicData_AlleleCounts_withmaf.csv` is not filtered for minor allele frequencies (MAF).
`FormattedFilteredGenomicData_AlleleCounts_withoutmaf.csv` is filtered for MAF < 1%.

**Dataset `FormattedFilteredGenomicData_AlleleFrequencies_withoutmaf.csv`**

Allele frequencies of the populations filtered for MAF < 1%. SNPs in columns and populations in rows.

  1.  `pop`: population.
  2. to 9818. SNPs ID.
  

&nbsp;

#### Imputed genomic data

Missing values were imputed based on the most common allele within the main gene pool of the genotype of concern (although we acknowledge that some genotypes had high admixture rates). See section 6 of the report `1_FormattingGenomicData.qmd`.

**Datasets `ImputedGenomicData_AlleleCounts_withmaf.csv` and `ImputedGenomicData_AlleleCounts_withmaf.csv`**

Imputed allele counts of each genotype (coded as 0, 1 or 2) . SNPs in rows and genotypes in columns.

  1. `snp_ID`: SNP ID.
  2. to 455. `ALT10` to `VER9`:  Clone (i.e., genotype) ID.
  
`ImputedGenomicData_AlleleCounts_withmaf.csv` is not filtered for MAF.
`ImputedGenomicData_AlleleCounts_withoutmaf.csv` is filtered for MAF < 1%.

**Datasets `ImputedGenomicData_AlleleFrequencies_withmaf.csv` and `ImputedGenomicData_AlleleFrequencies_withoutmaf.csv`**

Imputed allele frequencies of the populations. SNPs in columns and populations in rows.

  1.  `pop`: population.
  2. to 9818. SNPs ID.
  
`ImputedGenomicData_AlleleFrequencies_withmaf.csv` is not filtered for MAF.
`ImputedGenomicData_AlleleFrequencies_withoutmaf.csv` is filtered for MAF < 1%.

&nbsp;

#### Genomic position and other SNP IDs

**Datasets `ListSNPs_withmaf.csv` and `ListSNPs_withoutmaf.csv`**

Additional information on SNPs position on the genome and the different SNPs IDs used across studies/assays. SNPs are in rows.

  1. `snp_ID`: SNP ID used in the present study.
  2. `original_ID`: Original ID of the SNP.
  3. `affx_ID`: SNP ID in the Axiom assay.
  4. `infinium_ID`: SNP ID in the Infinium assay.
  5. `scaffold/contig`: contig on which the SNP is located. The term scaffold is also used as some SNPs were obtained from the alignment of NGS short-reads from a pseudoreference genome in *Pinus pinaster* that were called scaffolds (even though there are not really scaffolds).
  6. `genome_position`: position of the SNP on the scaffold/contig (in bp).
  7. `annotation`: SNP name (including possible alternative names).

`ListSNPs_withmaf.csv` includes SNPs with minor allele frequencies (MAF).
`ListSNPs_withoutmaf.csv` does not include SNPs with MAF < 1%.


&nbsp;

### Climatic data

Climatic data was extracted with the [Climate Downscaling Tool (ClimateDT)](https://www.ibbr.cnr.it/climate-dt/). The full name, units and calculation method of the climatic variables can be found on the ClimateDT website: https://www.ibbr.cnr.it/climate-dt/?action=fldlist. Note that the variable called mean summer precipitation (`MSP`) in ClimateDT was renamed summer precipitation (`SP`) in the present study. 

  - `tmn01` to `tmn12`: monthly minimum temperature (°C).
  - `tmx01` to `tmx12`: monthly maximum temperature (°C).
  - `prc01` to `prc12`: monthly total precipitation (mm).
  - `bio1`: mean annual temperature (°C). 
  - `bio2`: mean diurnal range (mean of monthly (maximum temperature - minimum temperature)) (°C).
  - `bio3`: isothermality (`bio2`/`bio7`) (×100) (index).
  - `bio4`: temperature seasonality (standard deviation ×100) (°C).
  - `bio5`: maximum temperature of the warmest month (°C).
  - `bio6`: minimum temperature of the coldest month (°C).
  - `bio7`: temperature annual range (`bio5`-`bio6`) (°C).
  - `bio8`: mean temperature of the wettest quarter (°C).
  - `bio9`: mean temperature of the driest quarter (°C).
  - `bio10`: mean temperature of the warmest quarter (°C).
  - `bio11`: mean temperature of the coldest quarter (°C).
  - `bio12`: annual precipitation (mm).
  - `bio13`: precipitation of the wettest month (mm).
  - `bio14`: precipitation of the driest month (mm).
  - `bio15`: precipitation seasonnality (coefficient of variation) (index).
  - `bio16`: precipitation of the wettest quarter (mm).
  - `bio17`: precipitation of the driest quarter (mm).
  - `bio18`: precipitation of the warmest quarter (mm).
  - `bio19`: precipitation of the coldest quarter (mm).
  - `Eref`: Hargreaves reference evaporation (mm).
  - `CMD`: Hargreaves climatic moisture deficit (mm).
  - `MCMT`:mean coldest month temperature (°C).
  - `SP`: summer (May to September) precipitation (mm).
  - `MWMT`: summer (May to September) precipitation (mm).
  - `GDD0`: degree-days above 0°C, chilling degree-days (°C * days).
  - `GDD5`: degree-days above 5°C, chilling degree-days (°C * days).
  - `GDD18`: degree-days above 18°C, chilling degree-days (°C * days).
  - `NFFD`: number of frost-free days (day).
  - `DMA`: De Martonne aridity index (index).
  - `EPQ`: Emberger pluviothermic quotient (index).
  - `RMT`: Rivas-Martinex thermic index (index).
  - `AHM`: annual heat moisture index (bio1+10)/(bio12/1000)) (°C/mm).
  - `SHM`: summer heat moisture index MWMT/(MSP/1000) (°C/mm).
  - `TD`: temperature difference MWMT-MCMT (°C).
  - `EMT`: extreme minimum temperature over 30 years (°C).
  - `bFFP`: Julian date on which FFP begins (day).
  - `eFFP`: Julian date on which FFP ends (day).
  - `FFP`: frost-free period (day).
  - `PAS`: precipitation as snow (mm).
  - `spi01` to `spi12`: monthly standardized precipitation index (index).
  - `spei01` to `spei12`: monthly standardized precipitation index (index).



**Dataset `ClimateDT_CommonGardens.csv`**

Annual climatic data (from 2010 to 2018) at the location of the five common gardens: Asturias (Spain), Bordeaux (France), Cáceres (Spain), Madrid (Spain) and Fundão (Portugal). Dataset generated in `2_CommonGardenData.qmd`.

  1. `cg`: common garden.
  2. `longitude` of the common garden.
  3. `latitude` of the common garden.
  4. `elevation` of the common garden.
  5. `year` between 2010 and 2018.
  6. to 104. climatic variables (see names above).

**Dataset `AnnualClimaticData_FromClimateDT_PopulationLocations.csv`**

Annual climatic data (from 1901 to 2021) at the location of the 34 CLONAPIN populations. Dataset generated in `0_FormattingPopulationCoordinatesElevationClimateData.qmd`.

  1. `pop`: population.
  2. `longitude` of the population.
  3. `latitude` of the population.
  4. `elevation` of the population.
  5. `year` between 1901 and 2021.
  6. to 103. climatic variables (see names above).
  
**Dataset `AveragedClimaticVariables_ReferencePeriod_PopulationLocations.csv`**

Averages of the climatic variables over the reference period 1901-1950 at the location of the 34 CLONAPIN populations. Dataset generated in `0_FormattingPopulationCoordinatesElevationClimateData.qmd`.

  1. `pop`: population.
  2. `longitude` of the population.
  3. `latitude` of the population.
  4. `elevation` of the population.
  5. to 64. climatic variables (see names above).

**Dataset `AveragedClimaticVariables_FuturePeriod_PopulationLocations.csv`**

Predicted averages of the climatic variables at the location of the 34 CLONAPIN populations over the period 2041-2070, under the scenario SSP3.7-0 and for five climatic general circulation models (GCMs): GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0 and UKESM1-0-LL (column `gcm`). Dataset generated in `0_FormattingPopulationCoordinatesElevationClimateData.qmd`.

  1. `pop`: population.
  2. `gcm`: general circulation model.
  3. to 26. climatic variables (see names above).
  
**Dataset `ClimateDT_NFIPlots_PastClimates.csv`**

Averaged climatic data  at the location of the NFI plots for the reference period 1901-1950. Dataset generated in `13_ValidationNFI.qmd`.


**Dataset `ClimateDT_NFIPlots_SurveyClimates.csv`**

Averaged climatic data  at the location of the NFI plots for the survey periods specific to each plot. Dataset generated in `13_ValidationNFI.qmd`.

&nbsp;

### National Forest Inventory (NFI) data

**Dataset `NFIdata_cleaned.csv`**

This dataset contains mortality data from the natural populations of the Spanish and French NFI plots in which maritime pines were recorded. See report `13_ValidationNFI.qmd`. Meaning of the columns:

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
  11. `mean_DBH`: mean DBH of the maritime pines in the NFI plot including dead trees, adults and saplings/seedlings (recruitment).
  12. `prop_dead`: proportion of dead maritime pines in the NFI plot.
  13. `annual_prop_dead`: annual proportion of maritime pine mortality in the NFI plot.

&nbsp;

### Common garden data

Datasets generated in `15_ValidationCommonGardens.qmd`.

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


***

## SCRIPTS


`0_FormattingPopulationCoordinatesElevationClimateData.qmd` 

  - Checking population information (coordinates and elevation data) from different sources (e.g. collected from different studies).
  - Extracting climatic data from [Climate Downscaling Tool (ClimateDT)](https://www.ibbr.cnr.it/climate-dt/) at the location of the populations.
  - Calculating the average of the climatic variables across the time periods of interest.

`1_FormattingGenomicData.Rmd` 

  - Formatting genomic data: converting letters (e.g. A/A, A/G) to numbers (0,1 or 2), and `---` to `NA`.
  - Filtering genomic data for monomorphic SNPs, minor allele counts (MAC), proportion of missing data per clone and per SNP, minor allele frequencies (MAF).
  - Estimating statistical correlations among SNPs and LD.
  - Determining SNPs position on the genome.
  - Exploring genomic data, e.g., number of SNPs/clones genotyped in each assay, Average and maximum number of missing valuesper clone.
  - Imputation of missing data.

`2_CommonGardenData.qmd`

  - Extracting climatic data from [ClimateDT](https://www.ibbr.cnr.it/climate-dt/) at the location of the common gardens .
  - Calculating the mean climate in each common garden between the planting date and the measurement date.

`3_CheckingPastFutureClimatesPopulationLocations.qmd` 

  - Comparing ClimateDT climatic data from point estimates (generated using scale-free downscaling) and extracted values from rasters.
  - Comparing the values of the climatic variables at the location of the populations under two different reference periods, i.e., 1901-1950 and 1961-1990.
  - Comparing the values of the climatic variables at the location of the populations under current and future climates (from five GCMs).
  
`4_ReduncancyAnalysis_VariancePartionning_IdentificationCandidateSNPs.Rmd`

  - Selection of the climatic variables based on their biological relevance for maritime pine, their contribution to the genetic variance using a RDA-based stepwise selection ([Capblancq and Forester 2021](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13722)) and the magnitude of their exposure to climate change.
  - Partitioning genomic variation among climate, neutral population genetic structure (accounted for with the main axes of a PCA) and geography (accounted for with population coordinates or distance-based Moran's eigenvector maps).
  - Identification of the outlier SNPs using Redundancy analysis (RDA); approach developed in [Capblancq et al. (2018)](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12906) and [Capblancq and Forester (2021)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13722)).

`5_GradientForest_IdentificationCandidateSNPs.qmd`

  - Identification of outlier SNPs with the Gradient Forest (GF) algorithm, using either raw allele frequencies (GF-raw) or allele frequencies after correction for population relatedness (GF-X), as described in [Fitzpatrick et al. (2021)](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13374) and [Capblancq et al. (2023)](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.18465).

`6_BaypassAnalysis_IdentificationCandidateSNPs.qmd` 

  - Identification of outlier SNPs with [BayPass](https://www1.montpellier.inra.fr/CBGP/software/baypass/index.html); approach described in [Gautier (2015)](https://academic.oup.com/genetics/article/201/4/1555/5930067?login=true)).

`7_LFMM_IdentificationCandidateSNPs_GenomicOffsetPredictions.qmd` 

  - Identification of outlier SNPs and estimation of the genetic gap (i.e. genomic offset) with the `LEA` R package ([Gain & Francois 2021](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13366)), which uses the latent factor mixed model (LFMM) approach ([Frichot et al. 2013](https://academic.oup.com/mbe/article/30/7/1687/972098?login=true); [Cayes et al. 2019](https://academic.oup.com/mbe/article/36/4/852/5290100?login=true)). Mathematical and theoretical details of the approach in [Gain et al. (2023)](https://www.biorxiv.org/content/10.1101/2023.01.02.522469v3).

`8_GeneratingSNPsets.qmd`

  - Identifying the common outlier SNPs across the different gene-environment association (GEA) methods.
  - Checking the genome position of the outlier SNPs; when some of them were located on the same scaffold/contig, only the SNP with the lower $p$-value in the RDA was kept in the final set of candidate SNPs.
  - Generating a set of control SNPs (with the same number of SNPs as in the set of candidate SNPs).

`9_GeneralizedDissimilarityModelling_GenomicOffsetPredictions.qmd`

  - Genomic offset predictions with the Generalized Dissimilarity Modelling (GDM) approach, as described in [Fitzpatrick & Keller (2015)](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12376), [Mokany et al. (2022)](https://onlinelibrary.wiley.com/doi/full/10.1111/geb.13459) and the [GDM website](https://mfitzpatrick.al.umces.edu/gdm/).

`10_GradientForest_GenomicOffsetPredictions.qmd` 

  - Genomic offset predictions with the Gradient Forest (GF) algorithm, approach described in [Fitzpatrick & Keller (2015)](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12376) and [Gougherty et al. (2021)](https://www.nature.com/articles/s41558-020-00968-6).
    
`11_RedundancyAnalysis_GenomicOffsetPredictions.qmd` 

  - Genomic offset predictions with Redundancy Analysis (RDA); approach described in [Capblancq and Forester (2021)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13722).

`12_ComparingGenomicOffsetPredictions.qmd` 

  - Comparing genomic offset predictions across the different methods (GF, GDM, LFMM and RDA), SNPs sets (control and candidate SNPs; and also all SNPs for LFMM) and GCMs.

`13_ValidationNFI.qmd` 

  - Filtering and exploring mortality data from the National Forest Inventory (NFI) plots of France and Spain; see [Changenet et al. (2021)](https://onlinelibrary.wiley.com/doi/abs/10.1111/geb.13301). 
  - Extracting climatic data at the location of the NFI plots. 
  - Calculating the average of the climatic variables for the reference period (1901-1950) and the inventory period (specific to each NFI plot).
  - Estimating the association between the genomic offset predictions and mortality rates in the NFI plots. 

`14_ValidationNFI_ModelComparison.qmd` 

  - Building and evaluating the accuracy of the Bayesian models used to estimate the association between genomic offset predictions and mortality rates in the NFI plots.

`15_ValidationCommonGardens.qmd` 

  - Calculating the climate transfer distances (CTDs), i.e., absolute climatic difference between the location of the populations and the common gardens.
  - Estimating the association between genomic offset predictions / CTDs and mortality and height data from the five clonal common gardens (CLONAPIN network).
  - Comparing the predictive ability of genomic offset predictions vs CTDs.
  
`16_PCAplotIDpops.qmd` 

  - Generating figures for the Supplementary Information based on the PCA of the control and candidate SNPs: screeplots and PCA plots with the ALT and ARM populations highlighted.

`RPackageCitations.qmd` 

  - Citations of the R packages used in the present study.

