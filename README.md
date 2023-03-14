
# Losses and gains of fallows impact farmland bird populations over three funding periods of the EU Common Agricultural Policy
[![DOI](https://zenodo.org/badge/410768453.svg)](https://doi.org/10.5281/zenodo.7643282)


This is the companion repository for the manuscript: Associations between farmland bird and fallow area at large scales: consistently positive over three CAP periods but moderated by landscape complexity

## Structure of the repo

* script/ in this subfolder are four R scripts:
	1. 01_model_fitting_rich.R: is the script that runs the species richness models
	2. 02_model_fitting_abund.R: is the script that runs the abundance models
	3. 03_ plotting_script.R: is the script that produces the figures shown in the main text of the manuscript
	4. XX_helper_function.R: a file with helper functions used in the other scripts

* data/ in this subfolder are different datasets to reproduce the analysis, metadata information are given in the README of the data subfolder:
	1. bird_rich.csv: dataset used for the analysis of bird richness together with the predictor variables
	2. bird_abund.csv: dataset used for the species-level analysis of bird abundance together with the predictor variables
	2. bird_info.csv: a dataset with species names and group info of the different bird species
	4. district_fallow.gpkg: district geometries with fallow data

## Questions / Issues

For questions or issues either open an issue directly on the repo or contact sebastian.klimek[at]thuenen.de

## License

The files in the data/ subfolder are shared under the [CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0/) license. The files in the script/ subfolder are shared under the [MIT](https://github.com/github/choosealicense.com/blob/gh-pages/LICENSE.md) license.