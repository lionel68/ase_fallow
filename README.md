
# Losses and gains of fallows impact farmland bird populations over three funding periods of the EU Common Agricultural Policy
[![DOI](https://zenodo.org/badge/410768453.svg)](https://zenodo.org/badge/latestdoi/410768453)


This is the companion repository for the manuscript: Associations between farmland bird and fallow area at large scales: consistently positive over three CAP periods but moderated by landscape complexity

## Structure of the repo

* script/ in this subfolder are two R scripts:
	1. 01_model_fitting_rich.R: is the script that runs the species richness models
	2. 02_model_fitting_abund.R: is the script that runs the abundance models
        3. 03_ plotting_script.R: is the script that produces the figures shown in the main text of the manuscript
        4. XX_helper_function.R: a file with helper functions used in the other scripts

* data/ in this subfolder are different datasets, metadata information are given in the README of the data subfolder:
	1. bird_data.csv: main dataset used for the analysis containing bird scaled abundance and richness together with the predictor variables
	2. bird_info.csv: a dataset with species names and group info of the different bird species
	3. district_fallow.csv: district-level data (fallow land area ...)
	4. district_shape.gpkg: district geometries

## Questions / Issues

For questions or issues either open an issue directly on the repo or contact sebastian.klimek[at]thuenen.de
