####################
# Script to pre-process the raw data
# output the data files used in the analysis
##################################################

# note that due to data sharing policy the raw data
# are NOT on the repo; this script is therefore
# just the documentation of the data pre-processing steps

# pre-processing and data generation of ASE analysis
setwd("~/PostDoc_Thunen/ASE-stuff/")
## load libraries
library(tidyverse)
library(sf)
library(raster) # check if needed


# bird data
bird <- read.csv("data/birddata/mhb_bird.csv")
bird_group <- read.csv("data/birddata/bird_fallow_group.csv", encoding = "UTF-8")
names(bird)[1] <- "routcode"
# change two bunting to edge-breeders
bird_group$group[grep("Emberiza",bird_group$species)] <- "edge"


# location of CBBS transects
mhb_xy <- st_read("data/geodata/MhB_plots_S2637.shp")
mhb_xy <-st_drop_geometry(mhb_xy)[,c("ROUTCODE","X_COORD", "Y_COORD")]
names(mhb_xy)[1] <- tolower(names(mhb_xy)[1])

# compute the abundance per groups
bird %>%
  filter(! Artname %in% c("Rotmilan", "Rauchschwalbe")) %>% # remove red kite and barn swallow
  group_by(Artname) %>%
  mutate(abundance_scale = scales::rescale(Abund, c(0, 100))) %>% # re-scale abundance between 0-100
  ungroup() %>%
  left_join(bird_group) %>%
  group_by(routcode, Jahr, group) %>%
  summarise(Abundance = round(sum(abundance_scale), 0)) %>% # compute summed scaled abundance per group
  pivot_wider(names_from = group,
              values_from = Abundance) -> bird_abu

# compute the abundance per species
bird %>%
  filter(! Artname %in% c("Rotmilan", "Rauchschwalbe")) %>% # remove red kite and barn swallow
  group_by(routcode, Jahr, Artname) %>%
  summarise(Abundance = sum(Abund)) -> bird_abu2

# compute species richness and shannon div
bird %>%
  filter(! Artname %in% c("Rotmilan", "Rauchschwalbe")) %>%
  group_by(routcode, Jahr) %>%
  mutate(p = Abund / sum(Abund)) %>%
  summarise(rich = sum(p[!is.na(p) & p != 0] ** 0),
            shan = exp(-sum(p * log(p), na.rm = TRUE))) %>%
  mutate(shan = ifelse(rich == 0, 0, shan)) -> bird_div

# put abundance and diversity together
bird_dd <- full_join(bird_abu, bird_div,
                     by = c("routcode", "Jahr"))

# bkr info (not used in the analysis)
bkr <- read.csv("data/landusedata/mhb_bkr.csv")


# edge length data 
edge <- read.csv("data/landusedata/mhb_1kmbuffer_edge_atkisswf.csv")
edge$edge_m[is.na(edge$edge_m)] <- 0 # replace NAs by 0
# put in m edge per ha
edge$edge_mha <- (edge$edge_m) / (3.14e6 * 1e-4)


# atkis data
atkis <- read.csv("data/landusedata/mhb_1kmbuffer_atkis.csv")
# reclassify the landcover
pro_atki <- function(data, col){
  
  #col2 <- enquo(col)
  data$type <- ""
  data$type[data$obb %in% c("Siedlung", "Verkehr")] <- "urban"
  data$type[data$obb == "Gewässer"] <- "water"
  data$type[data$obabez %in% c("Wald, Forst ", "Gehölz ")] <- "forest"
  data$type[data$obabez %in% c("Sumpf, Ried ", "Moor, Moos ")] <- "water"
  data$type[data$obabez %in% c("Ackerland ", "Grünland ")] <- "agriculture"
  data$type[data$obabez %in% c("Heide ", "Gartenland ", "Sonderkultur ", "Vegetationslose Fläche ")] <- "rest"
  
  data %>%
    group_by(routcode, type) %>%
    summarise(area_ha = sum({{ col }} * 1e-4)) -> out
  
  return(out)
}
atkis %>%
  pro_atki(atkis_sqm) -> atkis_dd

# make this in wide format for merging
atkis_dd %>%
  filter(type != "") %>%
  pivot_wider(names_from = type, 
              values_from = area_ha, 
              values_fill = 0) -> atkis_dd2

# put together
bird_dd %>%
  ungroup() %>%
  rename(year = Jahr) %>%
  filter(year %in% c(2007, 2010, 2016)) %>%
  left_join(mhb_xy, by = "routcode") %>%
  left_join(edge[,-1], by = "routcode") %>%
  left_join(atkis_dd2, by = "routcode") %>%
  filter(!is.na(urban)) -> bird_dd2

# add missing edge infos
bird_dd2$edge_m[is.na(bird_dd2$edge_m)] <- 0
bird_dd2$edge_mha[is.na(bird_dd2$edge_mha)] <- 0

# remove plots with less than 10% agricultural area (ATKIS)
bird_dd2 %>%
  mutate(prop_agri = agriculture / (agriculture + forest + rest + urban + water)) %>%
  filter(prop_agri >= 0.1) -> bird_dd3

# now get fallow land area per district from ASE
ase <- read.csv("data/landusedata/agraratlas_2003_2016.csv")
# the muni shapefile
muni <- st_read("data/geodata/gemeinde_correct.shp")
# get the district shapefile
district <- st_read("data/geodata/250_NUTS3.shp")
district <- st_transform(district, st_crs(muni))
# create the NUTS code for the district
muni$NUTS_code <- substr(muni$GUI_CODE, 1, 5)
muni <- arrange(muni, NUTS_code)
names(ase)[1] <- "GEM_CHAR_G"

# some muni have incorrect NUTS_code, save them
muni_miss <- muni[!(muni$NUTS_code %in% district$NUTS_CODE),]
ii <- st_intersection(muni_miss, district)
ii$area <- st_area(ii)
ii %>%
  group_by(GUI_CODE) %>%
  filter(area == max(area)) %>%
  arrange(NUTS_code) -> ii_dd
# replace in the dataset with the correct NUTS code
muni$NUTS_code[muni$GUI_CODE %in% ii_dd$GUI_CODE] <- as.character(ii_dd$NUTS_CODE)



# aggregate per district and compute percent fallow land
# dividing by arable area, or if arable area is null
# by utilized arable area
ase %>%
  filter(variable %in% c("SEFA", "ARAB", "GRAS")) %>%
  inner_join(st_drop_geometry(muni), by = "GEM_CHAR_G") %>%
  group_by(NUTS_code, allYEAR, variable) %>%
  summarise(value = sum(value, na.rm = TRUE) * 1000) %>%
  pivot_wider(names_from = variable, values_from = value) %>%
  mutate(fallow_per = (SEFA / (ARAB + GRAS)) * 100) %>%
  rename(year = allYEAR) %>%
  filter(year %in% c(2007, 2010, 2016) & !is.na(year)) %>%
  unite("NUTS_year", NUTS_code, year, remove = FALSE) -> ase_district

# get the intersection between Mhb plots and the districts
mhb_xy <- st_as_sf(mhb_xy, coords = 2:3)
st_crs(mhb_xy) <- st_crs(muni)

ii <- st_intersects(mhb_xy, district)
mhb_xy$NUTS_code <- district$NUTS_CODE[as.numeric(ii)]

# merge district-level infos with the bird data
bird_dd3 %>%
  inner_join(st_drop_geometry(mhb_xy), by = "routcode") %>%
  inner_join(ase_district, by = c("NUTS_code", "year")) %>%
  filter(fallow_per < 100) -> bird_dd4

# save this
write.csv(bird_dd4, "data/preprocessed/bird_fallow_v10.csv",row.names = FALSE)

# count how many bird infos are available per district and year
bird_dd4 %>%
  group_by(NUTS_code, NUTS_year) %>%
  summarise(n_bird = n()) %>%
  mutate(year = as.integer(substr(NUTS_year, 7, 10))) %>%
  full_join(ase_district[,-c(1)], by = c("NUTS_code", "year")) %>%
  mutate(fallow_per = ifelse(fallow_per > 100, 100, fallow_per),
         is_bird = ifelse(is.na(n_bird), 0, 1),
         NUTS_year = ifelse(is.na(NUTS_year), 
                            paste(NUTS_code, year, sep = "_"),
                            NUTS_year)) -> ase_dd


write.csv(ase_dd, "data/preprocessed/district_fallow.csv", row.names = FALSE)


# now for species-level abundance
# put together
bird_abu2 %>%
  rename(year = Jahr) %>%
  filter(year %in% c(2007, 2010, 2016)) %>%
  left_join(edge[,-1], by = "routcode") %>%
  left_join(atkis_dd2, by = "routcode") %>%
  #left_join(bkr[,-1]) %>%
  mutate(prop_agri = agriculture / (agriculture + forest + urban + water + rest)) %>%
  filter(prop_agri >= 0.1) %>%
  filter(!is.na(urban)) %>%
  inner_join(st_drop_geometry(mhb_xy), by = "routcode") %>%
  inner_join(ase_district, by = c("NUTS_code", "year")) -> bird_spe

# add missing edge infos
bird_spe$edge_m[is.na(bird_spe$edge_m)] <- 0

# put in english names
bird_name <- read.csv("data/birddata/bird_fallow_group.csv",
                      stringsAsFactors = FALSE, encoding = "UTF-8")

bird_spe %>%
  left_join(bird_name) %>%
  mutate(species = gsub(" ", "_", species)) -> bird_spe2

# save this
write.csv(bird_spe2, "data/preprocessed/bird_fallow_species.csv",row.names = FALSE)
write.csv(bird_spe2, "Z:/Dokumente/Bird-stuff/monvia_voegel/ro2_ase_fallow/data/preprocessed/bird_fallow_species.csv",
          row.names=TRUE)
