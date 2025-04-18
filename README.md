# Evaluating genomic offset predictions in a forest tree with high population genetic structure


Juliette Archambeau<sup>1,2</sup>, Marta Benito-Garzón<sup>1</sup>, Marina de-Miguel<sup>1,3</sup>, Alexandre Changenet<sup>1</sup>, Francesca Bagnoli<sup>4</sup>, Frédéric Barraquand<sup>5</sup>, Maurizio Marchi<sup>4</sup>, Giovanni G. Vendramin<sup>4</sup>, Stephen Cavers<sup>2</sup>, Annika Perry<sup>2</sup> and Santiago C. González-Martínez<sup>1</sup>

**1** INRAE, Univ. Bordeaux, BIOGECO, F-33610 Cestas, France

**2** UK Centre for Ecology \& Hydrology, Bush Estate, Penicuik, United Kingdom

**3** EGFV, Univ. Bordeaux, Bordeaux Sciences Agro, INRAE, ISVV, F-33882, Villenave d'Ornon, France

**4** Institute of Biosciences and BioResources, National Research Council, 50019 Sesto Fiorentino, Italy

**5** CNRS, Institute of Mathematics of Bordeaux, F-33400 Talence, France

***



This repository contains all the scripts that we played around with during our long journey to arrive at the final version of Archambeau *et al*. (2025). Let’s be honest—scientific papers often take us on many detours that lead absolutely nowhere! If you’re after a tidier dataset that includes only the code used in the final paper, then I’d suggest checking out the Dryad repository, which also contains the data used in this work: <https://datadryad.org/dataset/doi:10.5061/dryad.bnzs7h4jt>.


***

**Paper abstract**

Genomic offset models are increasingly popular tools for identifying populations at risk of maladaptation under climate change. These models estimate the extent of genetic change required for populations to remain adapted under future climate change scenarios, but face strong limitations and still lack broad empirical testing. Using 9,817 single nucleotide polymorphisms (SNPs) genotyped in 454 trees from 34 populations of maritime pine, a species with a marked population genetic structure, we found substantial variability across genomic offset predictions from different methods, SNP sets, and general circulation models. Using five common gardens, we mostly found positive associations between genomic offset predictions and mortality, as expected. However, contrary to our expectations, we observed very few negative monotonic associations between genomic offset predictions and height. Higher mortality rates were also observed in national forest inventory plots with high genomic offset, but only for some methods and SNP sets. The differing genomic offset patterns produced by the best-validated methods across the maritime pine range hindered drawing definitive conclusions for the species. Our study demonstrates the imperative of employing different methods and validating genomic offset predictions with independent data sources before using them as reliable metrics to inform conservation or management.

***

## REPORTS

The code (`.qmd` and `Rmd` files) used to generate the following `html` reports are in the folder `/reports`.

##### [0_FormattingPopulationCoordinatesElevationClimateData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/0_FormattingPopulationCoordinatesElevationClimateData.html)

  - Checking population information (coordinates and elevation data) from different sources (e.g. collected from different studies).
  - Extracting climatic data from [Climate Downscaling Tool (ClimateDT)](https://www.ibbr.cnr.it/climate-dt/) at the location of the populations.
  - Calculating the average of the climatic variables across the time periods of interest.
  - Attributing to each clone and population its main gene pool based on the ancestry coefficients estimated with the STRUCTURE software in [Jaramillo-Correa et al. (2015)](https://academic.oup.com/genetics/article/199/3/793/5935834).

##### [1_FormattingGenomicData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/1_FormattingGenomicData.html)

  - Formatting genomic data: converting letters (e.g. A/A, A/G) to numbers (0,1 or 2), and `---` to `NA`.
  - Filtering genomic data for monomorphic SNPs, minor allele counts (MAC), proportion of missing data per clone and per SNP, minor allele frequencies (MAF).
  - Estimating statistical correlations among SNPs and LD.
  - Looking at SNP positions on the genome (i.e., scaffold/contigs of the SNPs).
  - Exploring genomic data, e.g., number of SNPs/clones genotyped in each assay, average and maximum number of missing values per clone.
  - Imputation of missing data.

##### [2_CommonGardenData.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/2_CommonGardenData.html)

  - Extracting climatic data from [ClimateDT](https://www.ibbr.cnr.it/climate-dt/) at the location of the common gardens .
  - Calculating the mean climate in each common garden between the planting date and the measurement date.

##### [3_CheckingPastFutureClimatesPopulationLocations.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/3_CheckingPastFutureClimatesPopulationLocations.html)

  - Comparing ClimateDT climatic data from point estimates (generated using scale-free downscaling) and extracted values from rasters.
  - Comparing the values of the climatic variables at the location of the populations under two different reference periods, i.e., 1901-1950 and 1961-1990.
  - Comparing the values of the climatic variables at the location of the populations under current and future climates (from five general circulation models).
  

##### [4_ReduncancyAnalysis_VariancePartionning_IdentificationCandidateSNPs.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/4_ReduncancyAnalysis_VariancePartionning_IdentificationCandidateSNPs.html)

  - Selection of the climatic variables based on their biological relevance for maritime pine, their contribution to the genetic variance using a RDA-based stepwise selection ([Capblancq and Forester 2021](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13722)) and the magnitude of their exposure to climate change.
  - Partitioning genomic variation among climate, neutral population genetic structure (accounted for with the main axes of a PCA) and geography (accounted for with population coordinates or distance-based Moran's eigenvector maps).
  - Identification of the outlier SNPs using Redundancy analysis (RDA); approach developed in [Capblancq et al. (2018)](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.12906) and [Capblancq and Forester (2021)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13722).
  
Some figures generated in this report:
    
  *  <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/PCAplot.pdf" target="_blank">PCAplot.pdf</a> Principal component analysis performed on the population allele frequencies. The first three axes of the PCA are used to account for the population structure in the RDA analysis.
  *  <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/RDAsummary.pdf" target="_blank">RDAsummary.pdf</a> Summary statistics of the RDA models.
  *  <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/RDAplots.pdf" target="_blank">RDAplots.pdf</a> RDA plots colored by gene pool.
  *  <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/RDAplots_outliers_1.pdf" target="_blank">RDAplots_outliers_1.pdf</a> RDA plots with outliers following [Forester et al. (2018)](https://onlinelibrary.wiley.com/doi/full/10.1111/mec.14584?casa_token=IOrVgFSER0gAAAAA%3AsOlFDnBLnWtTdC-R6vi5pZiRwuzpP4GQyr8H9hVpVqxW0_3RXOV6bznLQx9deVCrYv80LokfqFvaGeY) (and the [associated vignette](https://popgen.nescent.org/2018-03-27_RDA_GEA.html)).
  *  <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/RDAplots_outliers_2.pdf" target="_blank">RDAplots_outliers_2.pdf</a> RDA plots with outliers and  Manhattan plots following [Capblancq and Forester (2021)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13722) (and the [associated Github repository](https://github.com/Capblancq/RDA-landscape-genomics)).

##### [5_GradientForest_IdentificationCandidateSNPs.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/5_GradientForest_IdentificationCandidateSNPs.html) 

  - Identification of outlier SNPs with the Gradient Forest (GF) algorithm, using either raw allele frequencies (*GF-raw*) or allele frequencies after correction for population relatedness (*GF-X*), as described in [Fitzpatrick et al. 2021](https://onlinelibrary.wiley.com/doi/10.1111/1755-0998.13374) and [Capblancq et al. 2023](https://nph.onlinelibrary.wiley.com/doi/full/10.1111/nph.18465). Note that only outlier SNPs identified with *GF-raw* were used to select the potential candidate SNPs for adaptation to climate, which were then used to calculate the genomic offset.


##### [6_BaypassAnalysis_IdentificationCandidateSNPs.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/6_BaypassAnalysis_IdentificationCandidateSNPs.html) 

  - Identification of outlier SNPs with [BayPass](https://www1.montpellier.inra.fr/CBGP/software/baypass/index.html); approach described in [Gautier (2015)](https://academic.oup.com/genetics/article/201/4/1555/5930067?login=true).


##### [7_LFMM_IdentificationCandidateSNPs_GenomicOffsetPredictions.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/7_LFMM_IdentificationCandidateSNPs_GenomicOffsetPredictions.html)

  - Identification of outlier SNPs and estimation of the genetic gap (i.e. genomic offset) with the `LEA` R package ([Gain & Francois 2021](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13366)), which uses the latent factor mixed model (LFMM) approach ([Frichot et al. 2013](https://academic.oup.com/mbe/article/30/7/1687/972098?login=true); [Cayes et al. 2019](https://academic.oup.com/mbe/article/36/4/852/5290100?login=true)). Mathematical and theoretical details of the approach in [Gain et al. (2023)](https://www.biorxiv.org/content/10.1101/2023.01.02.522469v3).


##### [8_GeneratingSNPsets.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/8_GeneratingSNPsets.html) 

  - Identifying common outlier SNPs across the different gene-environment association (GEA) methods.
  - Checking the genome position of the outlier SNPs; when some of them were located on the same scaffold/contig, only the SNP with the lower $p$-value in the RDA was kept in the final sets of candidate SNPs.
  - Mapping SNP position on the reference genome of *Pinus tabuliformis* [(Niu et al. 2022)](https://doi.org/10.1016/j.cell.2021.12.006).
  - Generating three sets of candidate SNPs and two sets of control SNPs (with the same number of SNPs as in the set with all candidate SNPs).
  

##### [9_GeneralizedDissimilarityModelling_GenomicOffsetPredictions.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/9_GeneralizedDissimilarityModelling_GenomicOffsetPredictions.html)

  - Genomic offset predictions with the Generalized Dissimilarity Modelling (GDM) approach, as described in [Fitzpatrick & Keller (2015)](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12376), [Mokany et al. (2022)](https://onlinelibrary.wiley.com/doi/full/10.1111/geb.13459) and the [GDM website](https://mfitzpatrick.al.umces.edu/gdm/).

##### [10_GradientForest_GenomicOffsetPredictions.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/10_GradientForest_GenomicOffsetPredictions.html)

  - Genomic offset predictions with the Gradient Forest (GF) algorithm, approach described in [Fitzpatrick & Keller 2015](https://onlinelibrary.wiley.com/doi/full/10.1111/ele.12376) and [Gougherty et al. 2021](https://www.nature.com/articles/s41558-020-00968-6).

  
##### [11_RedundancyAnalysis_GenomicOffsetPredictions.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/11_RedundancyAnalysis_GenomicOffsetPredictions.html)

  - Genomic offset predictions with Redundancy Analysis (RDA) and partial RDA; approach described in [Capblancq and Forester (2021)](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13722).


##### [12_ComparingGenomicOffsetPredictions.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/12_ComparingGenomicOffsetPredictions.html) 

  - Comparing genomic offset predictions under future climates across the different methods (GF, GDM, LFMM, RDA and pRDA), SNPs sets (the three sets of controls SNPs and the three sets of candidate SNPs) and the five general circulation models. The variability in genomic offset predictions, considering all combinations of SNP sets, methods, and GCMs, along with their relationships to climatic distances, can be visualized [in the following Shiny app](https://juliettearchambeau.shinyapps.io/GenomicOffsetPredictionVariabilityInMaritimePine/).

##### [13_ValidationNFI.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/13_ValidationNFI.html) 

  - Filtering and exploring mortality data from the National Forest Inventory (NFI) plots of France and Spain; see [Changenet et al. (2021)](https://onlinelibrary.wiley.com/doi/abs/10.1111/geb.13301). 
  - Extracting climatic data at the location of the NFI plots. 
  - Calculating the average of the climatic variables for the reference period (1901-1950) and the inventory period (specific to each NFI plot).
  - Estimating the association between the genomic offset predictions and mortality rates in the NFI plots (with regression coefficients or correlations).

##### [14_ValidationNFI_ModelComparison.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/14_ValidationNFI_ModelComparison.html) 

  - Building and evaluating the accuracy of the Bayesian models used to estimate the association between genomic offset predictions and mortality rates in the NFI plots.


##### [15_ValidationCommonGardens.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/15_ValidationCommonGardens.html) 

  - Calculating the climate transfer distances (CTDs), i.e., absolute climatic difference between the location of the populations and the common gardens.
  - Estimating the association between genomic offset predictions / CTDs and mortality and height data from the five clonal common gardens (CLONAPIN network).
  - Comparing the predictive ability of genomic offset predictions vs CTDs.
  - Correlations among GO predictions and climatic transfer distances in each common garden are shown here: <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/Corrplot_GOpreds_CTDs_asturias.pdf" target="_blank">Asturias (Spain)</a>, <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/Corrplot_GOpreds_CTDs_bordeaux.pdf" target="_blank">Bordeaux (France)</a>, <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/Corrplot_GOpreds_CTDs_caceres.pdf" target="_blank">Cáceres (Spain)</a>, <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/Corrplot_GOpreds_CTDs_madrid.pdf" target="_blank">Madrid (Spain)</a>, <a href="https://juliettearchambeau.github.io/GOPredEvalPinpin/Corrplot_GOpreds_CTDs_portugal.pdf" target="_blank">Fundão (Portugal)</a>.

##### [16_PCAplotIDpops.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/16_PCAplotIDpops.html)

  - Generating additional figures based on the PCA of the control and candidate SNPs: screeplots and PCA plots with the ALT and ARM populations highlighted.


## License

The code of this repository is under the [MIT license](https://github.com/JulietteArchambeau/GOPredEvalPinpin/blob/main/LICENSE.md)

### Disclaimer

The functions below come from other sources and may therefore only be reused under the licenses indicated by their authors:

  - `baypass_utils.R`. R functions associated with the *BayPass* software (Gautier 2015). Available at <https://forgemia.inra.fr/mathieu.gautier/baypass_public>. Under the license [CeCILL-B FREE SOFTWARE LICENSE AGREEMENT](https://forgemia.inra.fr/mathieu.gautier/baypass_public/-/blob/master/LICENSE).
  - the `window` function in the script `calc_avg_clim_var.R` comes from the [`dismo` R package](https://cran.r-project.org/web/packages/dismo/). Under the [GNU GENERAL PUBLIC LICENSE, Version 3, 29 June 2007](https://cran.r-project.org/web/licenses/GPL-3). 
  - `corpmat.R` from Statistical Tools for High-throughput Data Analysis & Visualization (STHDA) website: <http://www.sthda.com/english/wiki/visualize-correlation-matrix-using-correlogram>.
  - `detectoutliers.R` from Forester et al. (2018): <https://popgen.nescent.org/2018-03-27_RDA_GEA.html>.
  - `rdadapt.R` from Capblancq & Forester (2021): <https://github.com/Capblancq/RDA-landscape-genomics/blob/main/src/rdadapt.R>.
  - the functions `format_geno`, `run_gf_ind`, `extract_pvals` and `identify_GFoutliers` in `5_GradientForest_IdentificationCandidateSNPs.qmd` were inspired from Fitzpatrick et al. (2021) and the associated [github repository](https://github.com/fitzLab-AL/geneticOffsetR). 


## Software versions

Analyses were undertaken with R version 4.3.3. R package versions are shown at the end of each report and in [RPackageCitations.html](https://juliettearchambeau.github.io/GOPredEvalPinpin/RPackageCitations.html).

## References

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
