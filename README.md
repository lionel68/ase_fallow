# Losses and gains of fallows drive farmland bird populations over three funding periods of the EU Common Agricultural Policy
[![DOI](https://zenodo.org/badge/410768453.svg)](https://zenodo.org/badge/latestdoi/410768453)

This is the companion repository for the manuscript: Losses and gains of fallows drive farmland bird populations over three funding periods of the EU Common Agricultural Policy

## Structure of the repo

* script/ in this subfolder are two R scripts:
	1. 01_model_fitting.R: is the script that runs the models and performs the model checks
	2. 02_ plotting_script.R: is the script that produces the figures shown in the main text of the manuscript
	3. 03_model_fitting_species_lvl.R: the script to run the models at the level of the species, note that the raw data to run these models are not given in the repo
	4. 04_sensitivity_analysis.R: the script to drop the species one at a time and compare the estimated effect of fallow land. 

* data/ in this subfolder are different datasets, metadata information are given in the README of the data subfolder:
	1. bird_data.csv: main dataset used for the analysis containing bird scaled abundance and richness together with the predictor variables
	2. district_fallow.csv: district-level data (fallow land area ...)
	3. district_shape.gpkg: district geometries

## Questions / Issues

For questions or issues either open an issue directly on the repo or contact lionel.hertzog[at]thuenen.de
