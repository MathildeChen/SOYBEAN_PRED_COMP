# Comparison of methods to aggregate climate data to predict crop yield: an application to soybean

## Summary
This repository contains scripts supporting a paper aiming to identify the best data-driven approach to predict soybean yields from climate inputs and to examine the impact of climate data aggregation method on predictive performances. These analyses were first conducted at the global scale and then at the scale of a single countries, in order to examine the robustness of the conclusions. The country studied are the United States of America and Brazil, two major areas producing soybean. 

## Packages required: 

## Analyses step: 

### Step 0: Home-made functions used to derive the different climate predictors 

### Step 1: Data pre-processing and training dataset constitution
The dataset (not provided in this repository) is derived from the global dataset of historical yields for major crops (GDHY, Iisumi and Sakai, 2020), ERA5-land database (Hersbach et al., 2023), and SPAM dataset on irrigation practices (Yu et al., 2020) which are all freely availble (see references below). The data were on 3447 sites worldwide from 1981 to 2016. 

The scripts used to perform this step are: 
- **01_1_Yield_data.R**: selection of the sites with enought soybean production (>1% of the surface dedicated to soybean) and those with 0 soybean (located in regions with unfavorable conditions for crop production)
- **01_2_Climate_data.R**: derivation of daily, monthly, and seasonaly climate variables for: minimum and maximum temperatures; mean precipitation; solar radiation; vapor pressure deficit; evapotranspiration.
- **01_Dataset_preparation.R**: merging datasets together

References and acces to the dataset: 
- Historical yields: Iizumi, T., Sakai, T. The global dataset of historical yields for major crops 1981–2016. Sci Data 7, 97 (2020). https://doi.org/10.1038/s41597-020-0433-7 ; the dataset is accessible here: https://doi.pangaea.de/10.1594/PANGAEA.909132
- Historical climate data: Hersbach, H., Bell, B., Berrisford, P., Biavati, G., Horányi, A., Muñoz Sabater, J., Nicolas, J., Peubey, C., Radu, R., Rozum, I., Schepers, D., Simmons, A., Soci, C., Dee, D., Thépaut, J-N. (2023): ERA5 hourly data on single levels from 1940 to present. Copernicus Climate Change Service (C3S) Climate Data Store (CDS), DOI: 10.24381/cds.adbb2d47 (Accessed in February 2023) ; the data is accessible here: https://cds.climate.copernicus.eu/datasets/reanalysis-era5-single-levels?tab=download  
- Irrigation fraction: Yu, Q., You, L., Wood-Sichra, U., Ru, Y., Joglekar, A. K. B., Fritz, S., Xiong, W., Lu, M., Wu, W., and Yang, P.: A cultivated planet in 2010 – Part 2: The global gridded agricultural-production maps, Earth Syst. Sci. Data, 12, 3545–3572, https://doi.org/10.5194/essd-12-3545-2020, 2020. SPAM2010 can be downloaded via an open-data repository (DOI: https://doi.org/10.7910/DVN/PRFF8V; IFPRI, 2019).

## Step 2: Deriving the different climate predictors
The temporal resolutions derived are 
- _'Daily'_ = Cumulative daily values over the growing season were computed for each climate variable.
- _'Monthly'_ = Montly average of non-cumulative daily climate data
- _'Seasonal'_ = Average over the entire soybean growing season of non-cumulative daily climate data.
- _'Standardized seasonal'_ = Average of all rescaled variables (i.e., each climate variable had a mean of zero and a standard deviation of one).

The five different dimension-reduction techniques applied to cumulative daily data and monthly averages were:
- principal componant analysis (PCA),
- functional principal componant analysis (FPCA),
- multivariate principal componant analysis (MFPCA),
- partial least square regression (PLSR),
- functional partial least square regression (FPLSR)

See supplementary material of the publication for methodological and computational details. 

These steps are performed using the following scripts: 
- **02_Dimension_reduction.R**: uses the functions from the **00_Functions_dimensions_reduction.R** script to derive the different scores at the global scale, then in the US, then in Brazil; this step generally takes very long to run, so I parallelized the analyses to run different sets of computation (e.g., deriving all the FPCA scores, then all MFPCA scores etc.)
- **02_Dimension_reduction_check_recomputation.R**: script used to check that the scores were correctly computed; 
- **02_Dimension_reduction_merge_scores.R**: merge the scores produced from several round of analyses (due to analyses parallelization);


## Step 3: Models development and comparison
In this step we simultaneously evaluate: 

(i) the impact of the temporal resolution of climate data (i.e., daily, monthly, seasonal, or standardized seasonal), 

(ii) a large range of dimension reduction techniques to aggregate climate data (i.e., PCA, FPCA, MFPCA, PLSR, FPLSR, or no reduction) 

(iii) modeling techniques to predict soybean yields from climate data (multiple linear regression or random forest). 

Models' prediction accuracy were evaluated using two metrics were used (root mean square error and an equivalent of R², the model efficiency) through two cross-validation procedure unsuring good model transferability in time and in space. 

This is performed in the:
- **03_Models_prediction.R** script; 

## Step 4: Figures production
These scripts produce the figures and supplementary figures of the paper. 

## Authors: 

- Mathilde Chen (https://orcid.org/0000-0002-5982-2143) $1,2,3$ 

- Nicolas Guilpart (https://orcid.org/0000-0003-3804-0211), ${4}$  

- David Makowski (https://orcid.org/0000-0001-6385-3703) ${1}$


Affiliations:

$1$  Université Paris-Saclay, INRAE, AgroParisTech, UMR MIA PS, 91120 Palaiseau, France

$2$ CIRAD, UMR PHIM, F-34398 Montpellier, France

$3$ PHIM, Univ Montpellier, CIRAD, INRAE, Institut Agro, IRD, Montpellier, France

$4$ Université Paris-Saclay, AgroParisTech, INRAE, UMR Agronomie, 91120 Palaiseau, France

## Related publication
Published 3 May 2024 • © 2024 The Author(s). Published by IOP Publishing Ltd

Environmental Research Letters, Volume 19, Number 5

Citation: Mathilde Chen et al 2024 Environ. Res. Lett. 19 054049

DOI 10.1088/1748-9326/ad42b5
