# Data and code for the paper: 'Evaluating genomic offset predictions in a forest tree with high population genetic structure'

Juliette Archambeau<sup>1,2</sup>, Marta Benito-Garzón<sup>1</sup>, Marina de-Miguel<sup>1,3</sup>, Alexandre Changenet<sup>1</sup>, Francesca Bagnoli<sup>4</sup>, Frédéric Barraquand<sup>5</sup>, Maurizio Marchi<sup>4</sup>, Giovanni G. Vendramin<sup>4</sup>, Stephen Cavers<sup>2</sup>, Annika Perry<sup>2</sup> and Santiago C. González-Martínez<sup>1</sup>


**1** INRAE, Univ. Bordeaux, BIOGECO, F-33610 Cestas, France

**2** UK Centre for Ecology \& Hydrology, Bush Estate, Penicuik, United Kingdom

**3** EGFV, Univ. Bordeaux, Bordeaux Sciences Agro, INRAE, ISVV, F-33882, Villenave d'Ornon, France

**4** Institute of Biosciences and BioResources, National Research Council, 50019 Sesto Fiorentino, Italy

**5** CNRS, Institute of Mathematics of Bordeaux, F-33400 Talence, France

**Corresponding author:** Juliette Archambeau, juli.archambeau@gmail.com

***

**Paper abstract**

Genomic offset models are increasingly popular tools for identifying populations at risk of maladaptation under climate change. These models estimate the extent of genetic change required for populations to remain adapted under future climate change scenarios, but face strong limitations and still lack broad empirical testing. Using 9,817 single nucleotide polymorphisms (SNPs) genotyped in 454 trees from 34 populations of maritime pine, a species with a marked population genetic structure, we found substantial variability across genomic offset predictions from different methods, SNP sets, and general circulation models. Using five common gardens, we mostly found positive associations between genomic offset predictions and mortality, as expected. However, contrary to our expectations, we observed very few negative monotonic associations between genomic offset predictions and height. Higher mortality rates were also observed in national forest inventory plots with high genomic offset, but only for some methods and SNP sets. The differing genomic offset patterns produced by the best-validated methods across the maritime pine range hindered drawing definitive conclusions for the species. Our study demonstrates the imperative of employing different methods and validating genomic offset predictions with independent data sources before using them as reliable metrics to inform conservation or management.

***

**Author contributions**

SCG-M designed the experiment and supervised the curation of field data. MdM cleaned and formatted the phenotypic data in the common gardens. MM extracted the climatic values at the location of the populations and the National Forest Inventory plots, and provided the raster files to project the genomic offset predictions across the species range. F Bagnoli, SCG-M and GGV did the DNA extractions, production, cleaning and checking of the genomic data. AC cleaned and formatted mortality data from the National Forest Inventory plots. SCG-M, JA and MBG conceived the paper methodology. JA and F Barraquand built the Bayesian models and codes in the validation steps. JA conducted the data analyses and wrote the code available in the present repository.

***

# Data

## Genomic data

Fo genomic data filtering and formatting, see `2_FormattingGenomicData.qmd`.

### Population genetic structure

#### `PopulationStructureCorrea2015.csv`

Proportion of assignment to the six gene pools identified in Jaramillo-Correa et al. (2015). This dataset contains 523 genotypes but only 454 genotypes were used in the present study because the other genotypes had no genomic data or too much missing data. Meaning of the columns:

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

### Raw genomic data

#### `RawGenomicData.csv`

Genomic data before data filtering and formatting. NAs are indicated with `---`. Meaning of the columns:

  1. `clone`: clone ID.
  2. `assay`: assay in which the clone was genotyped, either the Infinium assay (`only_Inf`), the Axiom assay (`only_Affx`) or both assays (`both_Inf_Affx`).
  3. `snp_1` to `snp_14016`:  genotype for each of the 14,016 SNPs.


&nbsp;


### Filtered and formatted allele counts

#### `FormattedFilteredGenomicData_AlleleCounts_withmaf.csv` and `FormattedFilteredGenomicData_AlleleCounts_withoutmaf.csv`

`FormattedFilteredGenomicData_AlleleCounts_withmaf.csv` contain alleles not filtered for minor allele frequencies (MAF); and `FormattedFilteredGenomicData_AlleleCounts_withoutmaf.csv` containes alleles filtered for MAF < 1%.. 
These datasets are allele counts of each genotype (coded as 0, 1 or 2). SNPs in rows and genotypes in columns. Meaning of the columns:

  1. `snp_ID`: SNP ID.
  2. to 455. `ALT10` to `VER9`:  Clone (i.e., genotype) ID.

&nbsp;

### Filtered and formatted allele frequencies

#### `FormattedFilteredGenomicData_AlleleFrequencies_withoutmaf.csv`

Allele frequencies of the populations filtered for MAF < 1%. SNPs in columns and populations in rows. Meaning of the columns:

  1.  `pop`: population.
  2. to 9818. SNPs ID.
  

&nbsp;


### Imputed allele counts

Missing values were imputed based on the most common allele within the main gene pool of the genotype of concern (although we acknowledge that some genotypes had high admixture rates).

#### `ImputedGenomicData_AlleleCounts_withmaf.csv` and `ImputedGenomicData_AlleleCounts_withmaf.csv`

`ImputedGenomicData_AlleleCounts_withmaf.csv` is not filtered for MAF. `ImputedGenomicData_AlleleCounts_withoutmaf.csv` is filtered for MAF < 1%.

These datasets contain imputed allele counts for each genotype (coded as 0, 1 or 2) . SNPs in rows and genotypes in columns.Meaning of the columns:

  1. `snp_ID`: SNP ID.
  2. to 455. `ALT10` to `VER9`:  Clone (i.e., genotype) ID.
  

&nbsp;

### Imputed allele frequencies

#### `ImputedGenomicData_AlleleFrequencies_withmaf.csv` and `ImputedGenomicData_AlleleFrequencies_withoutmaf.csv`

`ImputedGenomicData_AlleleFrequencies_withmaf.csv` is not filtered for MAF. `ImputedGenomicData_AlleleFrequencies_withoutmaf.csv` is filtered for MAF < 1%.

Theses datasets contain imputed allele frequencies of the populations. SNPs in columns and populations in rows. Meaning of the columns:

  1.  `pop`: population.
  2. to 9818. SNPs ID.
  


&nbsp;

### SNP information

#### `SNPsInformation.csv`. 

Information of the different SNP IDs used across studies/assays, SNP position on the genome, SNP annotation and whether the SNPs have a minor allele frequency. SNPs are in rows. Meaning of the columns:

  1. `snp_ID`: SNP ID used in the present study.
  2. `original_ID`: Original ID of the SNP.
  3. `affx_ID`: SNP ID in the Axiom assay.
  4. `infinium_ID`: SNP ID in the Infinium assay.
  5. `scaffold/contig`: contig on which the SNP is located. The term scaffold is also used as some SNPs were obtained from the alignment of NGS short-reads from a pseudoreference genome in *Pinus pinaster* that were called scaffolds (even though there are not really scaffolds).
  6. `genome_position`: position of the SNP on the scaffold/contig (in bp).
  7. `annotation`: SNP name (including possible alternative names).
  8. `MAF_filtering`: `REMOVED` for SNPs with minor allele frequencies, which are removed for GEAs analyses but not for estimating the neutral population genetic structure; `KEPT` otherwise.

&nbsp;

#### `SNPsInformation_WithOutliers.csv`. 

This dataset is similar to `SNPsInformation.csv` (the first eight columns are identical) but includes additional information about the identification of SNPs by the different gene-environment association (GEA) methods (GF, RDA, pRDA, LFMM and BayPass), and about the inclusion of the SNPs in the final sets of control and candidate SNPs (see `7_GeneratingSNPsets.qmd` report). Meaning of the columns: 

  9. `RDA_outliers`: `TRUE` if the SNP was identified as an outlier in the RDA; `FALSE` otherwise.
  10. `pRDA_outliers`: `TRUE` if the SNP was identified as an outlier in the pRDA; `FALSE` otherwise.
  11. `GF_outliers`: `TRUE` if the SNP was identified as an outlier with the GF algorithm; `FALSE` otherwise.
  12. `LFMM_outliers`: `TRUE` if the SNP was identified as an outlier woth LFMM; `FALSE` otherwise.
  13. `BayPass_outliers`: `TRUE` if the SNP was identified as an outlier with BayPass; `FALSE` otherwise.
  14. `all_candidate_SNPs`: `TRUE` if the SNP was an outlier identified by at least one GEA method among RDA, pRDA, LFMM, BayPass and GF; `FALSE` otherwise.
  15. `common_candidate_SNPs`: `TRUE` if the SNP was an outlier identified by at least two GEA methods among RDA, pRDA, LFMM, BayPass and GF; `FALSE` otherwise.
  16. `candidate_SNPs_corrected_pop_structure`: `TRUE` if the SNP was an outlier identified by at least one GEA method correcting for population genetic structure (i.e., pRDA, LFMM or BayPass); `FALSE` otherwise.
  17. `random_control_SNPs`: `TRUE` if the SNP was a control SNP randomly sampled among the SNPs not identified by any GEA method; `FALSE` otherwise.
  18. `control_SNPs_matching_allele_freq`: `TRUE` if the SNP was a control SNP sampled among the SNPs not identified by any GEA method, and has a similar allele frequency to a candidate SNP; `FALSE` otherwise.
  19. `control_SNPs_without_NAs`: `TRUE` if the SNP had no mising data; `FALSE` otherwise.

&nbsp;

#### `maritime_pine_coord.txt`

SNP positions on the chromosome-scale reference genome of a related species (i.e., *Pinus tabuliformis*; Niu et al. 2022), as such a reference is not yet available for maritime pine. Meaning of the columns:

  1. `AxiomID`: SNP ID in the Axiom assay.
  2. `CHROM`: assigned chromosome of the SNP on the reference genome. Contig was indicated when no chromosome was assigned.
  3. `Absolute_POS`: position of the SNP on the reference genome.
  

&nbsp;
  

## Climatic data

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

This information on the climatic variables used in the present study is also available in the `ClimateDT_Metadata.csv` dataset. In this dataset, the first column contains the full names of the climatic variables, the second column lists their units, and the third column provides their codes.


&nbsp;


### Climate in the common gardens

#### `ClimateDT_CommonGardens.csv`

Annual climatic data (from 2010 to 2018) at the location of the five common gardens: Asturias (Spain), Bordeaux (France), Cáceres (Spain), Madrid (Spain) and Fundão (Portugal).

  1. `cg`: common garden.
  2. `longitude` of the common garden.
  3. `latitude` of the common garden.
  4. `elevation` of the common garden.
  5. `year` between 2010 and 2018.
  6. to 104. climatic variables (see names above).
  

&nbsp;

### Climate at the location of the populations

#### `AnnualClimaticData_FromClimateDT_PopulationLocations.csv`

Annual climatic data (from 1901 to 2021) at the location of the 34 CLONAPIN populations.

  1. `pop`: population.
  2. `longitude` of the population.
  3. `latitude` of the population.
  4. `elevation` of the population.
  5. `year` between 1901 and 2021.
  6. to 103. climatic variables (see names above).
  
&nbsp;

#### `AveragedClimaticVariables_ReferencePeriod_PopulationLocations.csv`

Averages of the climatic variables over the reference period 1901-1950 at the location of the 34 CLONAPIN populations.

  1. `pop`: population.
  2. `longitude` of the population.
  3. `latitude` of the population.
  4. `elevation` of the population.
  5. to 64. climatic variables (see names above).

&nbsp;

#### `AveragedClimaticVariables_FuturePeriod_PopulationLocations.csv`

Predicted averages of the climatic variables at the location of the 34 CLONAPIN populations over the period 2041-2070, under the scenario SSP3.7-0 and for five climatic general circulation models (GCMs): GFDL-ESM4, IPSL-CM6A-LR, MPI-ESM1-2-HR, MRI-ESM2-0 and UKESM1-0-LL (column `gcm`).

  1. `pop`: population.
  2. `gcm`: general circulation model.
  3. to 26. climatic variables (see names above).
  
&nbsp;


### Climate in the inventory plots

#### `ClimateDT_NFIPlots_PastClimates.csv`

Averaged climatic data  at the location of the NFI plots for the reference period 1901-1950.

  1. `plotcode`: code of the inventory plot.
  2. `longitude` of the plot.
  3. `latitude` of the plot.
  4. `elevation` of the plot.
  5. to 64. climatic variables (see names above).

&nbsp;

#### `ClimateDT_NFIPlots_SurveyClimates.csv`

Averaged climatic data  at the location of the NFI plots for the survey periods specific to each plot.

  1. `plotcode`: code of the inventory plot.
  2. `longitude` of the plot.
  3. `latitude` of the plot.
  4. `elevation` of the plot.
  5. to 64. climatic variables (see names above).

&nbsp;


## National Forest Inventory (NFI) data

#### `NFIdata_cleaned.csv`

This dataset contains mortality data from the natural populations of the Spanish and French NFI plots in which maritime pines were recorded. Meaning of the columns:

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

## Common garden data

#### `CommonGardendata_cleaned.csv`

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

&nbsp;

#### `HeightIntercepts_Archambeauetal2023.csv`

This dataset contains the population height intercepts calculated across the five common gardens in the model 1 of Archambeau et al. (2022). Meaning of the columns:

  1. `pop`: population ID.
  2. `height`: mean of the posterior distributions of the population varying intercepts.
  3. `std.error`: standard error of the mean of the posterior distributions of the population varying intercepts.
  4. and 5. `conf.low` and `conf.high`: credible intervals of the posterior distributions of the population varying intercepts.


***

# Results

## Genomic offset predictions and climatic distances at the location of the populations

#### `go_predictions_climdist.csv`

Dataset generated in `11_ComparingGenomicOffsetPredictions.qmd`. Meaning of the columns:

  1. `var_name`: name of the set of SNPs or climatic variable.
  2. `method_name`: name of the genomic offset method or *CD* for climatic distance.
  3. `var_code`: code of the set of SNPs or climatic variable.
  4. `pop`: population.
  5. `val`: value of the genomic offset or the climatic distance.
  6. `gcm`: general circulation model used for future climate (2041-2060; SSP 3.7-0).
  7. `main_gp_pop`: main gene pool of the population.
  8. `color_main_gene_pool`: color used in the figures for the main gene pool of each population.

&nbsp;

## Genomic offset predictions and climatic distances at the location of the inventory plots

#### `go_nfi_predictions.csv`

Dataset generated in `12_ValidationNFI.qmd`. Meaning of the columns:

  1. `snp_set_code`: code of the set of SNPs (or code of the climatic variables for the climatic distances).
  2. `snp_set_name`: name of the set of SNPs (or name of the climatic variables for the climatic distances).
  3. `plotcode`: code of the inventory plot.
  4. `method`: name of the genomic offset method or *CD* for climatic distances.
  5. `GO`: value of the genomic offset or the climatic distance.
  6. `method_snp_set_code`: code of the combination of the genomic offset method (or climatic distance) and the set of SNPs (or climatic variable).
  7. `method_snp_set_name`: full name of the combination of the genomic offset method (or climatic distance) and the set of SNPs (or climatic variable).

&nbsp;

## Genomic offset predictions and climatic transfer distances at the location of the common gardens

#### `go_cg_predictions.csv`

Dataset generated in `13_ValidationCommonGardens.qmd`. Meaning of the columns:

  1. `cg`: name of the common garden.
  2. `pop`: population.
  3. `input_code`: code of the set of SNPs (or code of the climatic variables for the climatic transfer distances).
  4. `varX`: value of the genomic offset or the climatic transfer distance.
  5. `method`: name of the genomic offset method or *CTD* for climatic transfer distances.
  6. `method_input_code`: code of the combination of the genomic offset method (or climatic transfer distance) and the set of SNPs (or climatic variable).
  7. `input_name`: name of the set of SNPs (or code of the climatic variables for the climatic transfer distances).
  8. `method_input_name`: full name of the combination of the genomic offset method (or climatic transfer distance) and the set of SNPs (or climatic variable).



***

# Scripts

## Main reports

Quarto (or rmarkdown) and html versions of these scripts are available in the DRYAD repository.

All data used in these scripts are available in the DRYAD repository, except for some large rasters used to project genomic offset predictions across the range of maritime pine. These climatic data are available on the ClimateDT website: https://www.ibbr.cnr.it/climate-dt/.

The objects computed in these scripts are stored in an `/outputs` folder that users will need to create on their own computers to run the scripts successfully.

`1_MainPopulationGenePool_CommonGardenClimates.qmd` 

  - Attributing to each clone and population its main gene pool based on the ancestry coefficients estimated with the STRUCTURE software in Jaramillo-Correa et al. (2015).
  - Calculating the average of the climatic variables between the planting and the measurement dates in the common gardens.

`2_FormattingGenomicData.Rmd` 

  - Formatting genomic data: converting letters (e.g. A/A, A/G) to numbers (0,1 or 2), and `---` to `NA`.
  - Filtering genomic data for monomorphic SNPs, minor allele counts (MAC), proportion of missing data per clone and per SNP, minor allele frequencies (MAF).
  - Determining SNPs position on the genome.
  - Exploring genomic data, e.g., number of SNPs/clones genotyped in each assay, Average and maximum number of missing values per clone.
  - Imputation of missing data.
  
`3_ReduncancyAnalysis_VariancePartionning_IdentificationCandidateSNPs.Rmd`

  - Selection of the climatic variables based on their biological relevance for maritime pine, their contribution to the genetic variance using a RDA-based stepwise selection (Capblancq and Forester 2021) and the magnitude of their exposure to climate change.
  - Partitioning genomic variation among climate, neutral population genetic structure (accounted for with the main axes of a PCA) and geography (accounted for with population coordinates or distance-based Moran's eigenvector maps).
  - Identification of the outlier SNPs using Redundancy analysis (RDA) and partial RDA; approach developed in Capblancq et al. (2018) and Capblancq and Forester (2021).

`4_GradientForest_IdentificationCandidateSNPs.qmd`

  - Identification of outlier SNPs with the Gradient Forest (GF) algorithm as described in Fitzpatrick et al. (2021) and Capblancq et al. (2023).

`5_BaypassAnalysis_IdentificationCandidateSNPs.qmd` 

  - Identification of outlier SNPs with [BayPass](https://www1.montpellier.inra.fr/CBGP/software/baypass/index.html); an approach described in Gautier (2015).

`6_LFMM_IdentificationCandidateSNPs_GenomicOffsetPredictions.qmd` 

  - Identification of outlier SNPs and estimation of the genetic gap (i.e., genomic offset) with the `LEA` R package (Gain & Francois 2021), which uses the latent factor mixed model (LFMM) approach (Frichot et al. 2013, Caye et al. 2019). Mathematical and theoretical details of the approach in Gain et al. (2023).

`7_GeneratingSNPsets.qmd`

  - Identifying common outlier SNPs across the different gene-environment association (GEA) methods.
  - Checking the genome position of the outlier SNPs; when some of them were located on the same scaffold/contig, only the SNP with the lower $p$-value in the RDA was kept in the final sets of candidate SNPs.
  - Mapping SNP position on the reference genome of *Pinus tabuliformis* (Niu et al. 2022).
  - Generating three sets of candidate SNPs and two sets of control SNPs (with the same number of SNPs as in the set with all candidate SNPs).

`8_RedundancyAnalysis_GenomicOffsetPredictions.qmd` 

  - Genomic offset predictions with Redundancy Analysis (RDA) and partial RDA; approach described in Capblancq and Forester (2021).
  - Projections of the genomic offset predictions across the maritime pine range. Note that the rasters used for these projections (i.e., one for each climatic variable, and each general circulation model, so 6 rasters for the 1901-1950 reference climate; and 6*5=30 rasters for the 2041-2060 future climates) are not available in the DRYAD repository but are available on the [ClimateDT website](https://www.ibbr.cnr.it/climate-dt/).
  
`9_GeneralizedDissimilarityModelling_GenomicOffsetPredictions.qmd`

  - Genomic offset predictions with the Generalized Dissimilarity Modelling (GDM) approach, as described in Fitzpatrick & Keller (2015), Mokany et al. (2022) and the [GDM website](https://mfitzpatrick.al.umces.edu/gdm/).

`10_GradientForest_GenomicOffsetPredictions.qmd` 

  - Genomic offset predictions with the Gradient Forest (GF) algorithm, approach described in Fitzpatrick & Keller (2015) and Gougherty et al. (2021).
    
`11_ComparingGenomicOffsetPredictions.qmd` 

  - Comparing genomic offset predictions across the different methods (GF, GDM, LFMM, RDA and pRDA), SNPs sets (control and candidate SNPs) and GCMs.

`12_ValidationNFI.qmd` 

  - Exploring mortality data from the National Forest Inventory (NFI) plots of France and Spain; see Changenet et al. (2021). 
  - Calculating the climatic distances between the reference climate (1901-1950) and the climate during inventories at the location of the inventory plots.
  - Estimating the association between the genomic offset predictions and mortality rates in the NFI plots (with regression coefficients or correlations). 

`13_ValidationCommonGardens.qmd` 

  - Calculating the climate transfer distances (CTDs), i.e., absolute climatic difference between the climates at the location of the populations and the common gardens.
  - Estimating the association between genomic offset predictions / CTDs and mortality and height data from the five clonal common gardens (CLONAPIN network).

`RPackageCitations.qmd` 

  - Citations of the R packages used in the present study.

## Additional functions

Additionnal functions used in the main scripts: 
  
  - `calc_avg_clim_var.R` - function to calculate average climatic variables over a specific period, based on the monthly climatic data from ClimateDT. In this function, the `window` function comes from the [`dismo` R package](https://cran.r-project.org/web/packages/dismo/). Under the [GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007](https://cran.r-project.org/web/licenses/GPL-3). 
  - `corpmat.R` - to compute a matrix of p-values; from Statistical Tools for High-throughput Data Analysis & Visualization (STHDA) website: <http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram>.
  - `detectoutliers.R` - to identify outliers based on their RDA loadings (from  Forester et al. 2018); from Forester et al. (2018): <https://popgen.nescent.org/2018-03-27_RDA_GEA.html>.
  - `kable_mydf.R` - my own function to print table using the `kableExtra` R package.
  - `make_eucli_plot.R` - my own function to plot genomic offset predicions against Euclidean climatic distances, using R base plots.
  - `make_go_map.R` - my own function to plot genomic offset predictions at specific locations.
  - `rdadapt.R` - to conduct a RDA based genome scan; from Capblancq & Forester (2021): <https://github.com/Capblancq/RDA-landscape-genomics/blob/main/src/rdadapt.R>.


## *Stan* code of the models

*Stan* code of the models used in the validation steps:

  - `ValidationNFI_stancode.stan` - model used to estimate the association between mortality rates and genomic offset predictions/climatic distances in the inventory plots.
  - `ValidationNFI_stancode_with_predictions.stan` - same model as `ValidationNFI_stancode.stan` but including predictions of mortality rates in the model outputs.
  - `ValidationNFI_stancode_with_predictions_withoutGO.stan` - same model as `ValidationNFI_stancode_with_predictions.stan` but without genomic offset predictions or climatic distances as predictors (model used to estimate correlations between the model residuals and genomic offset predictions/climatic distances).
  - `ValidationCommonGarden_BinomialMortalityModel.stan` - model used to estimate the association between mortality rates and genomic offset predictions/CTDs in the common gardens.
  - `ValidationCommonGarden_BinomialMortalityModel_WithoutGO.stan` - same model as `ValidationCommonGarden_BinomialMortalityModel.stan` but without genomic offset predictions or climatic distances as predictors (model used to estimate correlations between the model residuals and genomic offset predictions/climatic distances).
  - `ValidationCommonGarden_HeightModel.stan` - model used to estimate the association between height and genomic offset predictions/CTDs in the common gardens.
  - `ValidationCommonGarden_HeightModel_WithoutPredictor.stan` - same model as `ValidationCommonGarden_HeightModel.stan` but without genomic offset predictions or climatic distances as predictors, model used to estimate the proportion of variance explained in a model that does not include the predictor of interest.
  - `ValidationCommonGarden_HeightModel_WithoutPredictor_WithPredictions.stan` - same model as `ValidationCommonGarden_HeightModel_WithoutPredictor.stan` but including height predictions in the model outputs.

***

# License

The code of this repository is under the MIT license.

## MIT license

Copyright (c) 2024 Authors

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

***

# Software versions

Analyses were undertaken with R version 4.3.3. R package versions are shown at the end of each script and in `RPackageCitations.html`.

***


# References

Archambeau J, Benito Garzón M, Barraquand F, de Miguel M, Plomion C and González-Martı́nez SC (2022). Combining climatic and genomic data improves range-wide tree height growth prediction in a forest tree. The American Naturalist 200(4):E141–E159.

Capblancq T, Luu K, Blum MGB and Bazin E (2018). Evaluation of redundancy analysis to identify signatures of local adaptation. Molecular Ecology Resources 18(6):1223–1233.

Capblancq T and Forester BR (2021). Redundancy analysis: A Swiss army knife for landscape genomics. Methods in Ecology and Evolution 12(12):2298–2309.

Capblancq T, Lachmuth S, Fitzpatrick MC and Keller SR (2023). From common gardens to candidate genes: exploring local adaptation to climate in red spruce. New Phytologist 237(5):1590–1605.

Caye K, Jumentier B, Lepeule J and François O (2019). LFMM 2: Fast and accurate inference of gene-environment associations in genome-wide studies. Molecular Biology and Evolution 36(4):852–860.

Changenet A, Ruiz-Benito P, Ratcliffe S, Fréjaville T, Archambeau J, Porte AJ et al. (2021). Occurrence but not intensity of mortality rises towards the climatic trailing edge of tree species ranges in European forests. Global Ecology and Biogeography 30(7):1356–1374.

Fitzpatrick MC and Keller SR (2015). Ecological genomics meets community-level modelling of biodiversity: mapping the genomic landscape of current and future environmental adaptation. Ecology Letters 18(1):1–16.

Fitzpatrick MC, Chhatre VE, Soolanayakanahally RY and Keller SR (2021). Experimental support for genomic prediction of climate maladaptation using the machine learning approach Gradient Forests. Molecular Ecology Resources 21(8):2749–2765.

Forester BR, Lasky JR, Wagner HH, Urban DL (2018). Comparing methods for detecting multilocus adaptation with multivariate genotype-environment associations. Molecular Ecology 27(9):2215-2233.

Frichot E, Schoville SD, Bouchard G and François O (2013). Testing for associations between loci and environmental gradients using latent factor mixed models. Molecular Biology and Evolution 30(7):1687–1699.

Gain C and François O (2021). LEA 3: Factor models in population genetics and ecological genomics with R. Molecular Ecology Resources 21(8):2738–2748.

Gain C, Rhoné B, Cubry P, Salazar I, Forbes F, Vigouroux Y et al. (2023). A quantitative theory for genomic offset statistics. Molecular Biology and Evolution 40(6):msad140.

Gautier M (2015). Genome-wide scan for adaptive divergence and association with population-specific covariates. Genetics 201(4):1555–1579.

Gougherty AV, Keller SR and Fitzpatrick MC (2021). Maladaptation, migration and extirpation fuel climate change risk in a forest tree species. Nature Climate Change 11:166–171. 

Jaramillo-Correa JP, Rodrı́guez-Quilón I, Grivet D, Lepoittevin C, Sebastiani F, Heuertz M et al. (2015). Molecular proxies for climate maladaptation in a long-lived tree (*Pinus pinaster* Aiton, Pinaceae). Genetics 199(3):793–807.

Mokany K, Ware C, Woolley SN, Ferrier S and Fitzpatrick MC (2022). A working guide to harnessing generalized dissimilarity modelling for biodiversity analysis and conservation assessment. Global Ecology and Biogeography 31(4):802-821.

Niu S, Li J, Bo W, Yang W, Zuccolo A, Giacomello S, ... & Wu HX (2022). The Chinese pine genome and methylome unveil key features of conifer evolution. Cell, 185(1), 204-217.

***

# Version changes

**29-jan-2025.**

These updates were made in response to revisions following the first round of review in *The American Naturalist*.

The abstract was modified. 

The dataset `SNPsInformation_WithOutliers.csv` was also updated, here are the columns that were modified or added:

  14. `all_candidate_SNPs`: `TRUE` if the SNP was an outlier identified by at least one GEA method among RDA, pRDA, LFMM, BayPass and GF; `FALSE` otherwise.
  15. `common_candidate_SNPs`: `TRUE` if the SNP was an outlier identified by at least two GEA methods among RDA, pRDA, LFMM, BayPass and GF; `FALSE` otherwise.
  16. `candidate_SNPs_corrected_pop_structure`: `TRUE` if the SNP was an outlier identified by at least one GEA method correcting for population genetic structure (i.e., pRDA, LFMM or BayPass); `FALSE` otherwise.
  17. `random_control_SNPs`: `TRUE` if the SNP was a control SNP randomly sampled among the SNPs not identified by any GEA method; `FALSE` otherwise.
  18. `control_SNPs_matching_allele_freq`: `TRUE` if the SNP was a control SNP sampled among the SNPs not identified by any GEA method, and has a similar allele frequency to a candidate SNP; `FALSE` otherwise.

Reports and scripts were updated to estimate the genomic offset for six sets of SNPs: three sets of candidate SNPs and three sets of control SNPs. Additionally, the estimation of genomic offset using pRDA was incorporated. Following these new analyses, the scripts for generating tables and figures were also updated.


**10-apr-2025.**

These updates were made to format the data and code according to the standards of *The American Naturalist*.

Modifications of the datasets:

  - `SNPsInformation_WithOutliers.csv` was divided into two datasets: `SNPsInformation.csv`, which contains information only about SNP IDs, genomic positions, annotations, and MAF; and `SNPsInformation_WithOutliers.csv`, which is the same dataset but includes additional information about outlier identification and the inclusion of SNPs in the sets of candidate and control SNPs.
  - Three datasets containing the genomic offset predictions and climatic distances estimated in the present study were uploaded in the following datasets: `go_predictions_climdist.csv`, `go_nfi_predictions.csv` and `go_cg_predictions.csv`.
  - In the dataset `ClimateDT_NFIPlots_SurveyClimates.csv`, the plots have been reordered to match the order in the dataset `ClimateDT_NFIPlots_PastClimates.csv`.
  - The `ClimateDT_Metadata.csv` dataset was added to provide information on the climatic variables from ClimateDT.
  - The `maritime_pine_coord.txt` dataset was added to provide the SNP positions on a reference genome (i.e., the *Pinus tabuliformis* genome).

Modifications of the scripts:

  - All reports have been revised between versions V2 and V3 to ensure the code is reproducible and includes only the sections of code relevant to the paper.
  - `0_FormattingPopulationCoordinatesElevationClimateData.qmd` and `2_CommonGardenData.qmd` were merged into `1_MainPopulationGenePool_CommonGardenClimates.qmd`.
  - `1_FormattingGenomicData.Rmd` was renamed to `2_FormattingGenomicData.Rmd`.
  - `3_CheckingPastFutureClimatesPopulationLocations.qmd`, `14_ValidationNFI_ModelComparison.qmd` and `16_PCAplotIDpops.qmd` were removed because they did not contain code that was directly used in the paper.
  - `4_ReduncancyAnalysis_VariancePartionning_IdentificationCandidateSNPs.Rmd` was renamed to `3_ReduncancyAnalysis_VariancePartionning_IdentificationCandidateSNPs.Rmd`.
  - `5_GradientForest_IdentificationCandidateSNPs.qmd` was renamed to `4_GradientForest_IdentificationCandidateSNPs.qmd`.
  - `6_BaypassAnalysis_IdentificationCandidateSNPs.qmd` was renamed to `5_BaypassAnalysis_IdentificationCandidateSNPs.qmd`.
  - `7_LFMM_IdentificationCandidateSNPs_GenomicOffsetPredictions.qmd`  was renamed to `6_LFMM_IdentificationCandidateSNPs_GenomicOffsetPredictions.qmd`.
  - `8_GeneratingSNPsets.qmd` was renamed to `7_GeneratingSNPsets.qmd`.
  - `11_RedundancyAnalysis_GenomicOffsetPredictions.qmd` was renamed to `8_RedundancyAnalysis_GenomicOffsetPredictions.qmd`.
  - `12_ComparingGenomicOffsetPredictions.qmd` was renamed to `11_ComparingGenomicOffsetPredictions.qmd`.
  - `13_ValidationNFI.qmd`  was renamed to `12_ValidationNFI.qmd`.
  - `15_ValidationCommonGardens.qmd` was renamed to `13_ValidationCommonGardens.qmd`.
  - The names of the files containing the *Stan* code for the models have been updated for clarity, and the content of each file is now explicitly described in the README.
  - Some additional scripts were removed as they were not relevant for the paper: `make_high_go_pop_maps.R`, `project_adaptive_index.R`, `generate_scaled_nfi_clim_datasets.R`, `generate_scaled_clim_datasets.R`, `extract_clim_from_rasters.R`, `extract_climatedt_metadata.R` and `baypass_utils.R`.
