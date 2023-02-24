# ReadMe #

This repository contains supplemental data and codes for the manuscript "Stochastic in space and time: Feedback between bedrock river incision and orographic patterns in runoff variability" by Forte & Rossi submitted to the Journal of Geophysical Research - Earth Surface. Below are brief descriptions of categories of items included in this repository. Note that the majority of processed data and all individual model results are hosted in a separate [Zenodo Repository](https://doi.org/10.5281/zenodo.7665887).

## Repository Contents
### geospatial_codes
Processing codes for analyzing either WaterGAP3 or GAGES II datasets and/or derivatives of these datasets. Majority of results of these scripts are stored in the related [Zenodo Repository](https://doi.org/10.5281/zenodo.7665887).
### manuscript_figures
Codes for generating the main text and supplmental figures for the in review manuscript
### model_runs
Model run scripts for all of the presented models. Note that outputs of the model are provided in the related [Zenodo Repository](https://doi.org/10.5281/zenodo.7665887).
### stimpy
Python codes for the 1D stocthastic threshold incision model

## Notes
### Raw Data
This repository primarily includes processing scripts to generate the processed datasets or figures, but it does include some of the processed data in the form of various data tables that are used in some of the scripts. Larger derived processed data and model outputs are hosted in a separate repository available on the related [Zenodo Repository](https://doi.org/10.5281/zenodo.7665887). These codes make use of a variety of publically available datasets (e.g., WaterGAP3 from Water Resources Reanalysis v2, GAGES II, etc.) that we do not rehost as we do not have the rights to distribute these data (see following section). We do provide the processing scripts that produced the derived products and those derived products.
### Required Publically Hosted Data
To successfully reproduce the analyses presented here, a variety of publically available data is required. Generally, in each script that makes use of these data sets, they are stored in a central location (denoted as 'master_location' in the codes), which you should update to an appropriate location on your machine if you wish to rerun or modify any of these scripts. The specific data required are:
* gagesII timeseries and watershed metadata from the [USGS](https://cmerwebmap.cr.usgs.gov/catalog/item/5788f619e4b0d27deb389055)
* 30 second STRM+ global topography from [Wordclim v2.1](https://www.worldclim.org/data/worldclim21.html)
* global ET dataset derived from Worldclim v2.1 from [Figshare](https://figshare.com/articles/dataset/Global_Aridity_Index_and_Potential_Evapotranspiration_ET0_Climate_Database_v2/7504448/4)
* several global 15 second rasters from [Hydrosheds v1](https://www.hydrosheds.org/hydrosheds-core-downloads)
* various earth2observe global reanalysis products for 1980-1989 and 1990-1999 from the FTP server at [earth2observe](http://www.earth2observe.eu/)
### Languages and Dependencies
The majority of the repository is written in Python, but there are isolated processing scripts written in Matlab. Below we list dependencies for the respective sets of scripts.
#### Python Dependencies
* numpy
* scipy
* matplotlib
* pandas
* scikit-learn
* cartopy
* cmcrameri
* fiona
* geopandas
* rasterio
* netcdf4
* gdal
#### Matlab Dependencies
* TopoToolbox
* Topographic-Analysis-Kit
