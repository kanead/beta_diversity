# Shrew MEM analysis
# code comes from the following paper:
# Benone et al. (2020) How modified landscapes filter rare species and modulate
# the regional pool of ecological traits?

#-------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls()) # removes everything currently held in the R memory
graphics.off() # closes all open graphics windows
# to clear the screen click "control L"
# set.seed(123)

## Libraries:
# ************
library(tidyverse)
library(adespatial)
library(spdep)
library(adegraphics)
library(gridExtra)
library(vegan)
library(usdm) # for VIF analysis
library(ade4)

#-------------------------------------------------------------------------------------------------------------------------------------------
## Useful function:
# *****************
nb2ggplot <- function(nb, coord) {
  # 'coord' must be a matrix/dataframe with two columns (called "long" and "lat")
  # 'nb' is an object of class nb
  # take out the connections from the nb object and assign them the lat and long in a dataframe
  n <- length(attributes(nb$neighbours)$region.id)
  DA <- data.frame(
    from = rep(1:n, sapply(nb$neighbours, length)),
    to = unlist(nb$neighbours),
    weight = unlist(nb$weights)
  )
  DA <- cbind(DA, coord[DA$from, 1:2], coord[DA$to, 1:2])
  colnames(DA)[4:7] = c("long", "lat", "long_to", "lat_to")
  return(DA)
}
set.seed(2020)

#-------------------------------------------------------------------------------------------------------------------------------------------
#' load the locations
shrewsite <- read.table("data/shrews_sites1.csv", sep = ",", head = TRUE)
shrewlocs <- shrewsite %>% dplyr::select(Longitude, Latitude)
shrewlocs <- shrewlocs %>% rename(long = Longitude)
shrewlocs <- shrewlocs %>% rename(lat = Latitude)
shrew.xy <- shrewlocs %>% as.matrix

# drop 67 sites for functional beta diversity
shrewlocs.func <- shrewsite %>% filter(Location %ni% c("CMR_Kupe1", "ETH_Alatish", "ETH_Bale1", "ETH_Bale3", "ETH_Bale4", "ETH_Bale5", "ETH_Chebera", "ETH_Hugumburda",
                                                       "ETH_Hunkolo", "ETH_Jiren", "ETH_Kaka", "ETH_Simiens2", "ETH_Simiens3", "ETH_Simiens4", "ETH_Yerer", "KEN_Kasigau",
                                                       "KEN_Meru", "KEN_Mukogodo", "Ken_Nairobi1", "KEN_Shimba", "KEN_Turkana1", "KEN_Turkana2", "LBR_Nimba1", "MLI_Mopti",
                                                       "MOZ_Cabora", "MOZ_Gorongosa", "MRT_Chott", "MWI_Lengwe", "MWI_Liwonde", "MWI_Zomba", "NAM_Keetmanshoop", "NAM_Tumasberg",
                                                       "NGA_Gambari", "SEN_Sabodala", "SLE_Bumbuna1", "SLE_Bumbuna2", "SOM_Jubba", "SWZ_Edwaleni", "SWZ_Maguga", "SWZ_Mlawula",
                                                       "TCD_Zakouma", "TZA_Kilimanjaro3", "TZA_Kingu", "UGA_Rwenzori4", "ZAF_Asante", "ZAF_Baviaanskloof", "ZAF_Karkloof",
                                                       "ZAF_Karoo2", "ZAF_Karoo3", "ZAF_Karoo4", "ZAF_Karoo5", "ZAF_Karoo6", "ZAF_Knersvlakte", "ZAF_Mankwe", "ZAF_Mariepskop",
                                                       "ZAF_MtZebra", "ZAF_Nossob", "ZAF_Nyslvley", "ZAF_Pretoriuskop", "ZAF_Rolfontein", "ZAF_Satara", "ZAF_Sneeuberg",
                                                       "ZAF_Swartberg", "ZAF_Wakefield", "ZAF_WillemPretorius", "ZMB_Kafue1", "ZWE_Sengwa")) %>% 
  dplyr::select(Longitude, Latitude) %>% rename(long = Longitude) %>% rename(lat = Latitude)
shrew.xy.func <- as.matrix(shrewlocs.func) 

shrewsp1 <-
  read.table("data/Shrews_spp1.csv",
             sep = ",",
             head = TRUE,
             row.names = 1)
shrewsp <- ifelse(shrewsp1 > 0, 1, 0)  # convert to presence-absence

# get environmental (BIOCLIM) data
enviro.shrew1 <- read.csv("data/bioclim.shrew.tax.csv")
enviro.shrew <- enviro.shrew1[,3:23]

# Select only variables with low VIF (<2.2) based on allsites selection (see Bat_MEM_AM.R)
env.shrew.final <- enviro.shrew %>% dplyr::select(alt, alt_rough, bio_2, bio_3, 
                                              bio_15, bio_18, bio_19)

# bioclim for functional beta diversity sites (remove 67)
env.shrew.func <- enviro.shrew1[,3:24] %>% filter(Location %ni% c("CMR_Kupe1", "ETH_Alatish", "ETH_Bale1", "ETH_Bale3", "ETH_Bale4", "ETH_Bale5", "ETH_Chebera", "ETH_Hugumburda",
                                                              "ETH_Hunkolo", "ETH_Jiren", "ETH_Kaka", "ETH_Simiens2", "ETH_Simiens3", "ETH_Simiens4", "ETH_Yerer", "KEN_Kasigau",
                                                              "KEN_Meru", "KEN_Mukogodo", "Ken_Nairobi1", "KEN_Shimba", "KEN_Turkana1", "KEN_Turkana2", "LBR_Nimba1", "MLI_Mopti",
                                                              "MOZ_Cabora", "MOZ_Gorongosa", "MRT_Chott", "MWI_Lengwe", "MWI_Liwonde", "MWI_Zomba", "NAM_Keetmanshoop", "NAM_Tumasberg",
                                                              "NGA_Gambari", "SEN_Sabodala", "SLE_Bumbuna1", "SLE_Bumbuna2", "SOM_Jubba", "SWZ_Edwaleni", "SWZ_Maguga", "SWZ_Mlawula",
                                                              "TCD_Zakouma", "TZA_Kilimanjaro3", "TZA_Kingu", "UGA_Rwenzori4", "ZAF_Asante", "ZAF_Baviaanskloof", "ZAF_Karkloof",
                                                              "ZAF_Karoo2", "ZAF_Karoo3", "ZAF_Karoo4", "ZAF_Karoo5", "ZAF_Karoo6", "ZAF_Knersvlakte", "ZAF_Mankwe", "ZAF_Mariepskop",
                                                              "ZAF_MtZebra", "ZAF_Nossob", "ZAF_Nyslvley", "ZAF_Pretoriuskop", "ZAF_Rolfontein", "ZAF_Satara", "ZAF_Sneeuberg",
                                                              "ZAF_Swartberg", "ZAF_Wakefield", "ZAF_WillemPretorius", "ZMB_Kafue1", "ZWE_Sengwa")) %>% 
  dplyr::select(alt, alt_rough, bio_2, bio_3, bio_15, bio_18, bio_19)
0#-------------------------------------------------------------------------------------------------------------------------------------------
# call on some shapefiles for plotting maps

biogeog1 <- readOGR("data/Bioregions_Linder_clipped_simplified.shp")
#proj <- "+proj=laea +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
proj1 <- '+proj=longlat + datum=WGS84'
biogeog <- spTransform(biogeog1, CRS(proj1))
proj4string(biogeog) # check crs

#africa1 <- rnaturalearth::ne_countries(continent = "africa", returnclass = "sf")

#-------------------------------------------------------------------------------------------------------------------------------------------
# A few graph-based connectivity schemes (before removing any link):
# ******************************************************************
shrew.nbtri <- tri2nb(shrew.xy) # Delaunay triangulation
shrew.nbgab <- graph2nb(gabrielneigh(shrew.xy), sym = TRUE) # Gabriel graph
shrew.nbrel <- graph2nb(relativeneigh(shrew.xy), sym = TRUE) # Relative neighbourhood graph
shrew.nbmst <- mst.nb(dist(shrew.xy)) # minimum spanning tree

# Visualisation of the connections:
# ********************************
# The spdep package is designed to create object plotted with base R.
# Visualisation of the connection schemes using ggplot2.
# We first need to transform the nb object into listw objects:
shrew.nbtri_listw <- nb2listw(shrew.nbtri)
shrew.nbgab_listw <- nb2listw(shrew.nbgab)
shrew.nbrel_listw <- nb2listw(shrew.nbrel)
shrew.nbmst_listw <- nb2listw(shrew.nbmst)

shrew.DA_tri <- nb2ggplot(shrew.nbtri_listw, shrew.xy)
shrew.tri_g <- shrew.xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = shrew.DA_tri,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "darkred"
  ) +
  labs(title = "Delaunay triangulation") +
  theme_bw()

shrew.DA_gab <- nb2ggplot(shrew.nbgab_listw, shrew.xy)
shrew.gab_g <- shrew.xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = shrew.DA_gab,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "darkred"
  ) +
  labs(title = "Gabriel") +
  theme_bw()

shrew.DA_rel <- nb2ggplot(shrew.nbrel_listw, shrew.xy)
shrew.rel_g <- shrew.xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = shrew.DA_rel,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "chocolate"
  ) +
  labs(title = "Relative neighbourhood") +
  theme_bw()

shrew.DA_mst <- nb2ggplot(shrew.nbmst_listw, shrew.xy)
shrew.mst_g <- shrew.xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = shrew.DA_mst,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "goldenrod"
  ) +
  labs(title = "Min. span tree") +
  theme_bw()
grid.arrange(shrew.tri_g, shrew.gab_g, shrew.rel_g, shrew.mst_g, ncol = 2, nrow = 2)

# These commented out lines work in basic R but I didn't change the networks
# nbgab_edited <- edit.nb(nbgab, xy)
# nbrel_edited <- edit.nb(nbrel, xy)
# save(nbgab_edited, nbrel_edited, file = "edited_nb_objects.RData")

# Load the edited nb objects:
# ---------------------------
# load("edited_nb_objects_shrews.RData")
# # Visualise the edited Gabriel graph and relative neighbourhood graphs:
# nbgab_listw_edited <- nb2listw(nbgab_edited)
# nbrel_listw_edited <- nb2listw(nbrel_edited)
# 
# 
# # 'nbdists' provides Euclidean distances along the links in a list of the same form as the nb list
# distgab <- nbdists(nbgab_edited, xy)
# # We want to replace these distance by the fluvial Euclidean distances, to better
# # represent the biology of the studied organisms:
# # The neighbours of each site are in ith element of the nb object list.
# # For example:
# nbgab_edited[[2]] # site 2 is connected to sites 86, 98

# Binary forms (no weights added to the connections):
shrew.nbtri_edited_b <- nb2listw(shrew.nbtri) # changed names here
shrew.nbgab_edited_b <- nb2listw(shrew.nbgab) # changed names here
shrew.nbrel_edited_b <- nb2listw(shrew.nbrel) # changed names here
shrew.nbmst_edited_b <- nb2listw(shrew.nbmst) # changed names here

#-------------------------------------------------------------------------------------------------------------------------------------------
# Construction of the list of candidate spatial weighting matrices:
# *****************************************************************
# Two candidates:
shrew.candidates <- list(shrew.tri_b = shrew.nbtri_edited_b,
                         shrew.gab_b = shrew.nbgab_edited_b,
                         shrew.rel_b = shrew.nbrel_edited_b,
                         shrew.nbm_b = shrew.nbmst_edited_b)

shrew.select <-
  listw.select(
    shrewsp,
    shrew.candidates,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )

# now save this object 'shrew.select' so that we don't have to rerun this!
save(shrew.select, file = "data/shrew.select.RData")
load("data/shrew.select.RData")

# Optimised selected SWM:
shrew.select$best.id

# Summarises the candidate SWMs:
shrew.select$candidates

# Summarises MEM var. selected and Moran's I and p-val in residuals after inclusion of the MEM var.
shrew.select$best$summary

# MEM variables to use in models and simulations:
# ***********************************************
shrew.MEM.sel <- shrew.select$best$MEM.select

# detrend latitude and longitude first
shrew.mods <- rda(shrewsp ~ long + lat, data = as.data.frame(shrew.xy), trace = FALSE)
shrew.res <- residuals(shrew.mods)

# now run this thru the same analysis on trended data (from above)
shrew.select.res <-
  listw.select(
    shrew.res,
    shrew.candidates,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )
shrew.MEM.sel.res <- shrew.select.res$best$MEM.select
shrew.select.res$candidates
save(shrew.select.res, file = "data/shrew.select.res.RData")
load("data/shrew.select.res.RData")

# now run this thru the same analysis but with dudi.pco of beta diversity (shrew.beta.pco = taxonomic Bsim)
# beta diversity pco calculated in 'aBetadiversity MEMs_Final.R'
shrew.select.beta.tax <-
  listw.select(
    shrew.beta.pco$li,
    shrew.candidates,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )
shrew.MEM.sel.beta.tax <- shrew.select.beta.tax$best$MEM.select
shrew.select.beta.tax$candidates
save(shrew.MEM.sel.beta.tax, file = "data/shrew.MEM.sel.beta.tax.RData")
load("data/shrew.MEM.sel.beta.tax.RData")

# now run this thru the same analysis but with dudi.pco of beta diversity (shrew.func.beta.pco = functional Bsim)
# beta diversity pco calculated in 'aBetadiversity MEMs_Final.R'
# but first need a new network because functional diversity dropped two sites

shrew.nbtri.func <- nb2listw(tri2nb(shrew.xy.func)) # Delaunay triangulation
shrew.nbgab.func <- nb2listw(graph2nb(gabrielneigh(shrew.xy.func), sym = TRUE)) # Gabriel
shrew.nbrel.func <- nb2listw(graph2nb(relativeneigh(shrew.xy.func), sym = TRUE)) # Relative neighbourhood
shrew.nbmst.func <- nb2listw(mst.nb(dist(shrew.xy.func))) # minimum spanning tree

# Four candidates
shrew.candidates.func <- list(shrew.tri_b = shrew.nbtri.func,
                            shrew.gab_b = shrew.nbgab.func,
                            shrew.rel_b = shrew.nbrel.func,
                            shrew.nbm_b = shrew.nbmst.func)

shrew.select.beta.func <-
  listw.select(
    shrew.beta.func.pco$li,
    shrew.candidates.func,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )
shrew.MEM.sel.beta.func <- shrew.select.beta.func$best$MEM.select
shrew.select.beta.func$candidates
save(shrew.MEM.sel.beta.func, file = "data/shrew.MEM.sel.beta.func.RData")
load("data/shrew.MEM.sel.beta.func.RData")

# now run this thru the same analysis but with dudi.pco of beta diversity (shrew.beta.phyl.pco = phylogenetic Bsim)
# beta diversity pco calculated in 'aBetadiversity MEMs_Final.R'
shrew.select.beta.phyl <-
  listw.select(
    shrew.beta.phyl.pco$li,
    shrew.candidates,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )
shrew.MEM.sel.beta.phyl <- shrew.select.beta.phyl$best$MEM.select
shrew.select.beta.phyl$candidates
save(shrew.MEM.sel.beta.phyl, file = "data/shrew.MEM.sel.beta.phyl.RData")
load("data/shrew.MEM.sel.beta.phyl.RData")

#-------------------------------------------------------------------------------------------------------------------------------------------
# Visualisation of the selected MEM variables:
# ********************************************

# s.value() not working, but I've created the MEMs in ggplot (see below)
s.value(
  shrew.xy,
  shrew.MEM.sel[, c(1:ncol(shrew.MEM.sel))],
  ppoint.cex = 0.6,
  symbol = "circle",
  xlim = c(min(shrew.xy[, 1]) - 0.3, max(shrew.xy[, 1]) + 0.3),
  ylim = c(min(shrew.xy[, 2]) - 0.3, max(shrew.xy[, 2]) + 0.3)
)

# view in ggplot (taken from 'Moran Eigenvector Map spatial structure.R')
shrew.mem.long1 <- shrew.MEM.sel %>% mutate('Site' = as.factor(seq(1, nrow(shrew.MEM.sel)))) %>% as.data.frame()
shrew.mem.long2 <- shrew.mem.long1 %>% gather(MEM, Value, MEM5:MEM29)
shrew.mem.long2$Longitude <- shrewsite$Longitude
shrew.mem.long2$Latitude <- shrewsite$Latitude

shrew.mem.long <- shrew.mem.long2 %>% group_by(Longitude, Latitude, MEM) %>% mutate('Mean' = mean(Value)) %>% 
  filter(MEM %in% c("MEM5","MEM10","MEM7","MEM13","MEM11","MEM3","MEM6","MEM8","MEM4","MEM2","MEM15",
                    "MEM9","MEM26","MEM16","MEM12","MEM31","MEM42","MEM37","MEM19","MEM1","MEM70",
                    "MEM44","MEM68","MEM23","MEM28","MEM25","MEM20","MEM40","MEM60","MEM57","MEM51",
                    "MEM22","MEM27","MEM21","MEM29"))

shrew.mem.long$MEM = factor(shrew.mem.long$MEM, levels=c("MEM5","MEM10","MEM7","MEM13","MEM11","MEM3","MEM6","MEM8","MEM4","MEM2","MEM15",
                                                         "MEM9","MEM26","MEM16","MEM12","MEM31","MEM42","MEM37","MEM19","MEM1","MEM70",
                                                         "MEM44","MEM68","MEM23","MEM28","MEM25","MEM20","MEM40","MEM60","MEM57","MEM51",
                                                         "MEM22","MEM27","MEM21","MEM29"))

shrew.mem.fort <- fortify(shrew.mem.long)
biogeog.fort <- merge(fortify(biogeog), as.data.frame(biogeog), by.x="id", by.y=0) # this is needed to add "Region"
#biogeog.fort <- merge(fortify(africa1), as.data.frame(africa1), by.x="sovereignt", by.y=0) # this is needed to add "Region"

shrew.mem.graph <- ggplot() +
  # geom_polygon(data=biogeog.fort, aes(x=long, y = lat, group=group), colour='black', 
  #              fill = c('white')) + # I'm failing to give each region a different colour - 'fill = piece' not working
  #coord_fixed(xlim = c(-72,-45), ylim=c(35,70), ratio=1.2)+
  facet_wrap(~MEM, ncol=11, nrow=5) +
  #facet_grid(~MEM) +
  #facet_wrap_paginate(~MEM, page=2) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.title = element_blank())+
  geom_point(data=shrew.mem.fort, aes(x = Longitude, y = Latitude, size=Value, fill=Value), shape=22) +
  theme_bw()+theme(legend.position = "right",
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(y= 'Latitude') +  
  labs(x= 'Longitude') +
  theme(axis.text.x=element_text(colour="white")) +
  theme(axis.text.y=element_text(colour="white")) +
  scale_fill_gradient2(low='black', mid='white', high="red", midpoint = 0) #+ # Adam, how can I manually set the size of the dots based on the value of 'Value'? I want to exaggerate the size difference

#-------------------------------------------------------------------------------------------------------------------------------------------
# Visualisation of the species composition multiscale spatial patterns:
# *********************************************************************
rda <- rda(shrewsp ~ . , data = MEM.sel)
# Permutation test (although the correct p-value is the adjusted one obtained with listw.select)
anova.cca(rda)

# Permutation test of each constraint axis (RDA axes):
test_by_axis <- anova.cca(rda, by = "axis")
nb_ax <- length(which(test_by_axis$"Pr(>F)" <= 0.05))
paste("The first", nb_ax, "RDA axes are significant.", sep = " ")

# Adjusted R-squared of the model:
RsquareAdj(rda)$adj.r.squared # must correspond to AdjR2 obtained with select$best$summary

sc <- scores(
  rda,
  choices = c(1:nb_ax),
  display = "lc",
  scaling = 1
)
s.value(
  xy,
  sc[, c(1:nb_ax)],
  ppoint.cex = 0.6,
  symbol = "circle",
  xlim = c(min(xy[, 1]) - 0.3, max(xy[, 1]) + 0.3),
  ylim = c(min(xy[, 2]) - 0.3, max(xy[, 2]) + 0.3)
)

#-------------------------------------------------------------------------------------------------------------------------------------------
# now do variance partitioning using varpart() from vegan
# what I'm not sure about is which axis to use from the rda. I've used 'wa' under 'CCA'
# shrew.vp <- varpart(shrew.rda$CCA$wa, enviro.shrew, shrew.MEM.sel) # can use < transfo = "hel" > 
# # but gives identical result, which isn't surprising (Hellinger transformation is for abundance data)
# plot(shrew.vp, bg = c(3, 5), Xnames = c("environment", "spatial"))

# repeat but using species-sites data (ie 'batsp') - this should be correct now
shrew.vp1 <- varpart(shrewsp, env.shrew.final, shrew.MEM.sel) 
plot(shrew.vp1, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))

# now repeat with detrended lat and long (see above) - should lat-long be here? Seems useless
shrew.vp1.res <- varpart(shrew.res, env.shrew.final, shrew.MEM.sel.res)
plot(shrew.vp1.res, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))

# now repeat with pco of beta diversity (tax =taxonomic, func = functional)
shrew.vp1.beta.tax <- varpart(shrew.beta.tax.pco$li, env.shrew.final, shrew.MEM.sel.beta.tax)
shrew.vp1.beta.func <- varpart(shrew.beta.func.pco$li, env.shrew.func, shrew.MEM.sel.beta.func)
shrew.vp1.beta.phyl <- varpart(shrew.beta.phyl.pco$li, env.shrew.final, shrew.MEM.sel.beta.phyl)

par(mar = c(1.5, 3.5, 1.5, 3.5)) # c(bottom, left, top, right)
par(mfrow = c(3, 1))
plot(shrew.vp1.beta.tax, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))
plot(shrew.vp1.beta.func, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))
plot(shrew.vp1.beta.phyl, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))
par(mfrow = c(1, 1))
par(mar = c(4, 4, 1.5, 0.1)) 

#-------------------------------------------------------------------------------------------------------------------------------------------
# now divide into broad, meso, and fine scale and analyse relative contributions
# of environment and space separately at each scale

# Now repeat but this time use the species sites data (original analysis deleted from above)
# broad scale = MEM: 5, 10, 7, 13 (four)
# meso scale = MEM: 15, 9, 16, 19 (four)
# low scale = MEM: 70, 44, 40, 51 (four)

shrew.vp.broad <- varpart(shrewsp, env.shrew.final, shrew.MEM.sel[, c(1,2,3,4)])  
shrew.vp.meso <- varpart(shrewsp, env.shrew.final, shrew.MEM.sel[, c(11,12,14,19)])  
shrew.vp.fine <- varpart(shrewsp, env.shrew.final, shrew.MEM.sel[, c(21,22,28,31)]) 

par(mar = c(1, 1.5, 1.5, 1)) 
par(mfrow = c(3, 1))
plot(shrew.vp.broad, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))
plot(shrew.vp.meso, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))
plot(shrew.vp.fine, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))
par(mfrow = c(1, 1))
par(mar = c(4, 4, 1.5, 0.1)) 
