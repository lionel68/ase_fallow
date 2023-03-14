# Metadata information

Below are metadata information on the files present in the data/ subfolder.


* bird_rich.csv
	* NUTS_anon: anonymized code for the district (NUTS level 3)
	* rich: farmland bird species richness (of the 24 focal species)
	* edge_std: z-standardized edge density between woody features (forest, hedges ...) and agricultural land (arable plus grassland) (unitless)
	* agri_std_: z-standardized area (unitless) covered by agricultural land (arable plus grassland)
	* fallow_std: z-standardized square root of the percent of fallow area on agricultural area (arable plus grassland)
	* year: the year of sampling of the bird and agricultural data (fallow and total agricultural area). Note that edge density is assumed constant across the years.

* bird_abund.csv
	* NUTS_anon: anonymized code for the district (NUTS level 3)
	* Abundance: abundance proxy from the CBBS monitoring
	* edge_std: z-standardized edge density between woody features (forest, hedges ...) and agricultural land (arable plus grassland) (unitless)
	* agri_std_: z-standardized area (unitless) covered by agricultural land (arable plus grassland)
	* fallow_std: z-standardized square root of the percent of fallow area on agricultural area (arable plus grassland)
	* year: the year of sampling of the bird and agricultural data (fallow and total agricultural area). Note that edge density is assumed constant across the years.
	* Species.name: english name of the bird species

* district_fallow.gpkg
	* year: the year of sampling from the agricultural statisticx
	* year_cat*: re-naming of the year variable
	* fallow_cut: proportion of fallow area over the total agricultural area in classes
	* fallow_d_cut: changes in proportion of fallow area between the sampling years in classes
	* geom: the geometry

