# ReadMe #

This is the ReadMe for the GitHub repository that contains supplemental data and codes for the in review manuscript "Stochastic in space and time: Feedback between bedrock river incision and orographic patterns in runoff variability" by Forte & Rossi submitted to the Journal of Geophysical Research - Earth Surface. Below are brief descriptions of categories of items included in this repository.

## Repository Contents
### geospatial_codes
Processing codes for analyzing either WaterGAP3 or GAGES II datasets and/or derivatives of these datasets.
### manuscript_figures
Codes for generating the main text and supplmental figures for the in review manuscript
### model_runs
Model run scripts for all of the presented models. Note that outputs of the model are provided in the related [Figshare Repository].
### stimpy
Python codes for the 1D stocthastic threshold incision model

## Notes
### Raw Data
This repository primarily includes processing scripts to generate the processed datasets or figures, but it does include some of the processed data in the form of various data tables that are used in some of the scripts. Larger derived processed data and model outputs are hosted in a separate repository available on the related [Figshare Repository]. These codes make use of a variety of publically available datasets (e.g., WaterGAP3 from Water Resources Reanalysis v2, GAGES II, etc.) that we do not rehost as we do not have the rights to distribute these data. We do provide the processing scripts that produced the derived products and those derived products.
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