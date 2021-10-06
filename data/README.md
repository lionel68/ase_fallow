# Metadata information

* bird_data.csv
	* NUTS_code: code for the district (NUTS level 3), note that these codes do not match with the other dataset on purpose
	* NUTS_year: NUTS_code together with the year of measurement
	* edge: scaled abundance for farmland bird species breeding in the edges of fallows
	* field: scaled abundance for farmland bird species breeding in fallow fields
	* low: scaled abundance for farmland bird species not breeding in or at the edge of fallows
	* rich: farmland bird species richness
	* edge_mha: edge density between woody features (forest, hedges ...) and agricultural land (arable plus grassland) in meter per hectaes
	* agriculture: area (ha) covered by agricultural land (arable plus grassland)
	* fallow_per: percent of fallow area on agricultural area (arable plus grassland)

* district_fallow.csv
	* NUTS_code: official NUTS code for the district, this code correspond to the code in the geometry file
	* NUTS_year: NUTS_code together with the year of measurement
	* year: year of the measurement
	* ARAB: arable land area in the district in ha
	* GRAS: grassland area in the district in ha
	* SEFA: fallow land area in the district in ha
	* fallow_per: percent of fallow land in the district computed as: SEFA / (GRAS + ARAB)

* district_shape.gpkg, geometries of the district with NUTS_code