# pre-processing and data generation of ASE analysis
setwd("~/PostDoc_Thunen/ASE-stuff/")
## load libraries
library(tidyverse)
library(sf)


# bird data
bird <- read.csv("data/birddata/mhb_bird.csv")
bird_group <- read.csv("data/birddata/bird_fallow_group.csv")
names(bird)[1] <- "routcode"

# location of CBBS transects
mhb_xy <- st_read("data/geodata/MhB_plots_S2637.shp")
mhb_xy <-st_drop_geometry(mhb_xy)[,c("ROUTCODE","X_COORD", "Y_COORD")]
names(mhb_xy)[1] <- tolower(names(mhb_xy)[1])

# compute the abundance per groups, [add weighing max abu]
bird %>%
  filter(! Artname %in% c("Rotmilan", "Rauchschwalbe")) %>% # remove red kite and barn swallow
  group_by(Artname) %>%
  mutate(abundance_scale = round(scales::rescale(Abund, c(0, 100)), 0)) %>%
  ungroup() %>%
  left_join(bird_group) %>%
  group_by(routcode, Jahr, group) %>%
  summarise(Abundance = sum(abundance_scale)) %>%
  pivot_wider(names_from = group,
              values_from = Abundance) -> bird_abu

# compute the abundance per species
bird %>%
  filter(! Artname %in% c("Rotmilan", "Rauchschwalbe")) %>% # remove red kite and barn swallow
  group_by(routcode, Jahr, Artname) %>%
  summarise(Abundance = sum(Abund)) -> bird_abu2

# compute species richness and shnnon div
bird %>%
  filter(! Artname %in% c("Rotmilan", "Rauchschwalbe")) %>%
  group_by(routcode, Jahr) %>%
  mutate(p = Abund / sum(Abund)) %>%
  summarise(rich = sum(p[!is.na(p) & p != 0] ** 0),
            shan = exp(-sum(p * log(p), na.rm = TRUE))) %>%
  mutate(shan = ifelse(rich == 0, 0, shan)) -> bird_div

bird_dd <- full_join(bird_abu, bird_div,
                     by = c("routcode", "Jahr"))

# fallow data
fallow <- read.csv("data/landusedata/mhb_1kmbuffer_fallowatkis.csv")

# bkr data
bkr <- read.csv("data/landusedata/mhb_bkr.csv")

# swf data 
#swf <- read.csv("data/landusedata/mhb_1kmbuffer_swf.csv")
# put together LWF and PWF
# swf %>%
#   group_by(routcode) %>%
#   summarise(swf_sqm = sum(swf_sqm)) -> swf_dd
# # edge length data
# edge <- read.csv("data/landusedata/mhb_1kmbuffer_edge.csv")
# edge length data combining ATKIS and SWF
edge <- read.csv("data/landusedata/mhb_1kmbuffer_edge_atkisswf.csv")
edge$edge_m[is.na(edge$edge_m)] <- 0 # replace NAs by 0

# atkis data
atkis <- read.csv("data/landusedata/mhb_1kmbuffer_atkis.csv")
pro_atki <- function(data, col){
  
  #col2 <- enquo(col)
  data$type <- ""
  data$type[data$obb %in% c("Siedlung", "Verkehr")] <- "urban"
  data$type[data$obb == "Gewässer"] <- "water"
  data$type[data$obabez %in% c("Wald, Forst ", "Gehölz ")] <- "forest"
  data$type[data$obabez %in% c("Sumpf, Ried ", "Moor, Moos ")] <- "water"
  data$type[data$obabez %in% c("Gartenland ", "Sonderkultur ", "Ackerland ")] <- "agriculture"
  data$type[data$obabez %in% c("Grünland ", "Heide ")] <- "grassland"
  
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

# load the intensity data
intensity <- read.csv("data/landusedata/intensity_costs.csv", sep = ";", dec =",")
names(intensity)[1] <- "id"
# remove 1% highest value for DAIR, GRAN, OGRL
rem_99 <- function(x,p=0.99){
  q99 <- quantile(x, probs = p)
  x_out <- ifelse(x>q99, q99, x)
  return(x_out)
}
# extend per usage
intensity %>%
  pivot_wider(id_cols = c(id, ha),
              names_from = agg_verfahren,
              values_from = kosten_pro_ha) %>%
  mutate(across(3:7, ~ifelse(is.na(.x), 0, .x))) %>% # fill na rows with 0
  mutate(across(5:7, ~rem_99(.x))) %>% # the 1% top values in DAIR, OGRL and GRAN are replaced
  pivot_longer(3:7, names_to="type",
               values_to = "cost_per_ha") %>%
  mutate(cost = cost_per_ha * ha) %>% # re-compute cost
  group_by(id) %>%
  summarise(ges_kost_ha = sum(cost)) -> int_dd

# this is sent to the db for the estimation based on ATKIS agricultural land
write.csv(int_dd, "data/landusedata/intensity_total_cost.csv", row.names = FALSE)

# load back
intensity <- read.csv("data/landusedata/mhb_intensity_atkis.csv")

# put together
bird_dd %>%
  rename(year = Jahr) %>%
  filter(year %in% c(2007, 2010, 2016)) %>%
  left_join(mhb_xy, by = "routcode") %>%
  left_join(intensity, by = "routcode") %>%
  left_join(edge[,-1], by = "routcode") %>%
  left_join(atkis_dd2, by = "routcode") %>%
  left_join(fallow, by = c("routcode", "year")) %>%
  left_join(bkr[,-1]) %>%
  filter(!is.na(urban)) %>%
  filter(!is.na(fallow_ha)) -> bird_dd2

# add missing bkr infos
bird_dd2$bkr[bird_dd2$routcode == "ni130"] <- unique(bird_dd2$bkr[bird_dd2$routcode == "ni129"])
bird_dd2$bkr[bird_dd2$routcode == "sh42"] <- unique(bird_dd2$bkr[bird_dd2$routcode == "sh43"])
# add missing edge infos
bird_dd2$edge_m[is.na(bird_dd2$edge_m)] <- 0

# remove plots with less than 10% agricultural area (ATKIS)
bird_dd2 %>%
  mutate(prop_agri = agriculture / (agriculture + forest + grassland + urban + water)) %>%
  mutate(fallow_ha = fallow_ha * 1000) %>%
  filter(prop_agri >= 0.1) -> bird_dd3

# save this
write.csv(bird_dd3, "data/preprocessed/bird_fallow_v5.csv",row.names = FALSE)

# now for species-level abundance
# put together
bird_abu2 %>%
  rename(year = Jahr) %>%
  filter(year %in% c(2007, 2010, 2016)) %>%
  left_join(intensity, by = "routcode") %>%
  left_join(edge[,-1], by = "routcode") %>%
  left_join(atkis_dd2, by = "routcode") %>%
  left_join(fallow, by = c("routcode", "year")) %>%
  left_join(bkr[,-1]) %>%
  filter(!is.na(urban)) %>%
  filter(!is.na(fallow_ha)) -> bird_spe

# remove plots with less than 10% agricultural area (ATKIS)
bird_spe %>%
  mutate(prop_agri = agriculture / (agriculture + forest + grassland + urban + water)) %>%
  mutate(fallow_ha = fallow_ha * 1000) %>%
  filter(prop_agri >= 0.1) -> bird_spe2

# add missing bkr infos
bird_spe2$bkr[bird_spe2$routcode == "ni130"] <- unique(bird_spe2$bkr[bird_spe2$routcode == "ni129"])
bird_spe2$bkr[bird_spe2$routcode == "sh42"] <- unique(bird_spe2$bkr[bird_spe2$routcode == "sh43"])
# add missing edge infos
bird_spe2$edge_m[is.na(bird_spe2$edge_m)] <- 0

# put in english names
bird_name <- read.csv("data/birddata/bird_fallow_group.csv",
                      stringsAsFactors = FALSE)

bird_spe2 %>%
  left_join(bird_name) %>%
  mutate(species = gsub(" ", "_", species)) -> bird_spe2

# save this
write.csv(bird_spe2, "data/preprocessed/bird_fallow_species.csv",row.names = FALSE)
write.csv(bird_spe2, "Z:/Dokumente/Bird-stuff/monvia_voegel/ro2_ase_fallow/data/preprocessed/bird_fallow_species.csv",
          row.names=TRUE)
