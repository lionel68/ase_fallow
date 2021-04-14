# additional graphs building up

# analysis script
setwd("~/PostDoc_Thunen/ASE-stuff/")

# load libraries
library(tidyverse)
library(sf)


# load data
bird <- read.csv("data/preprocessed/bird_fallow.csv")
bird$land <- substr(bird$routcode, 1, 2)

mhb <- st_read("data/geodata/MhB_plots_S2637.shp")
names(mhb)[17] <- "routcode"
mhb <- st_drop_geometry(mhb)

swf <- read.csv("data/landusedata/mhb_1kmbuffer_swf.csv")
swf$percent <- swf$sum / (pi*1000**2)

de <- st_read("data/geodata/DE_mainland.shp")

mhb %>%
  mutate(land = substr(routcode, 1, 2)) %>%
  filter(!(land %in% c("be", "hb", "hh") & is.na(BUND))) %>%
  select(routcode, lat, long_) %>%
  inner_join(swf) %>%
  group_by(code) %>%
  mutate(swf_cat = Hmisc::cut2(percent, g=6))-> mhb_swf

gg_1 <- ggplot() +
  geom_sf(data=de, fill = NA) +
  geom_point(data=subset(mhb_swf, code==1), aes(x=long_, y=lat, color=swf_cat)) +
  scale_color_brewer(palette = "BuGn",
                     labels = c("0 - 0.5", "0.5 - 1.2", "1.2 - 2.0", "2.0 - 2.9",
                                "2.9 - 4.4", "4.4 - 16.4"),
                     name = "Percent Linear\nWoody Features") +
  labs(x="", y="") +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())

gg_2 <- ggplot() +
  geom_sf(data=de, fill = NA) +
  geom_point(data=subset(mhb_swf, code==2), aes(x=long_, y=lat, color=swf_cat)) +
  scale_color_brewer(palette = "BuGn",
                     labels = c("0 - 0.3", "0.3 - 0.7", "0.7 - 1.2", "1.2 - 1.8",
                                "1.8 - 2.7", "2.7 - 10.1"),
                     name = "Per mil Patchy\nWoody Features") +
  labs(x="", y="") +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())

gg_3 <- ggplot() +
  geom_sf(data=de, fill = NA) +
  geom_point(data=subset(mhb_swf, code==3), aes(x=long_, y=lat, color=swf_cat)) +
  scale_color_brewer(palette = "BuGn",
                     labels = c("0 - 0.3", "0.3 - 0.7", "0.7 - 1.2", "1.2 - 2.0",
                                "2.0 - 3.3", "3.3 - 20.5"),
                     name = "Percent Additional\nWoody Features") +
  labs(x="", y="") +
  theme(axis.text = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.background = element_blank())  

gg_a <- gridExtra::grid.arrange(gg_1, gg_2, gg_3, ncol=3)
ggsave("figures/swf_mhb.png", gg_a, width = 8, height = 6)

atkis <- read.csv("data/landusedata/mhb_1kmbuffer_atkis.csv")

atkis$type <- ""
atkis$type[atkis$obb %in% c("Siedlung", "Verkehr")] <- "urban"
atkis$type[atkis$obb == "Gewässer"] <- "water"
atkis$type[atkis$obabez %in% c("Wald, Forst ", "Gehölz ")] <- "forest"
atkis$type[atkis$obabez %in% c("Sumpf, Ried ", "Moor, Moos ")] <- "water"
atkis$type[atkis$obabez %in% c("Gartenland ", "Sonderkultur ", "Ackerland ")] <- "agriculture"
atkis$type[atkis$obabez %in% c("Grünland ", "Heide ")] <- "grassland"

atkis %>%
  filter(type != "") %>%
  group_by(routcode, type) %>%
  summarise(percent = sum(sum) / (pi * 1000 ** 2)) %>%
  pivot_wider(names_from = type,
              values_from = percent,
              values_fill = 0) -> atkis_per

# try some clustering approaches
dd <- dist(atkis_per[,-1])

cc <- hclust(dd, method="ward.D")
gr <- cutree(cc, k = 6)
atkis_per$group <- factor(gr,
                          labels = c("agriculture", "mixed agriculture /\nforest",
                                     "mixed agriculture /\ngrassland / wet areas",
                                     "urban", "forest", "mixed forest/\ngrassland"))


atkis_per %>%
  pivot_longer(2:6) -> pp

gg_box <- ggplot(pp, aes(x=factor(group), y=value)) +
  geom_boxplot() +
  facet_wrap(~name) +
  labs(x="Clustering groups",
       y="Percent land-use")
# group 1: agriculture (arable crop)
# group 2: mixed agriculture / forest
# group 3: mixed agriculture / grassland / wet areas
# group 4: urban areas
# group 5: forest
# group 6: mixed forest / grassland


# plot this in space
mhb %>%
  mutate(land = substr(routcode, 1, 2)) %>%
  filter(!(land %in% c("be", "hb", "hh") & is.na(BUND))) %>%
  select(routcode, lat, long_) %>%
  inner_join(atkis_per) %>%
  mutate(group_f = factor(group,
                          labels = c("agriculture", "mixed agriculture /\nforest",
                                     "mixed agriculture /\ngrassland / wet areas",
                                     "urban", "forest", "mixed forest/\ngrassland"))) -> aa

gg_map <- ggplot() +
  geom_sf(data=de, fill = NA) +
  geom_point(data=aa, aes(x=long_, y=lat, color=group_f)) +
  scale_color_manual(values=c("lightgoldenrod1",
                              "lightgoldenrod4",
                              "aquamarine3",
                              "darkorange3",
                              "darkgreen",
                              "darkolivegreen3"),
                     name = "Clustering groups") +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.line = element_blank(),
        panel.background = element_blank()) +
  labs(x="", y="")

gg_atkis <- gridExtra::grid.arrange(gg_map, gg_box, ncol=2)  
ggsave("figures/atkis_clustering.png", gg_atkis, width = 16, height = 10)
