# Rat MEM analysis
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
library(rgdal)
library(sf)  # to project lat/long for variogram
#library(ggforce) # for facet_wrap_paginate() but didn't need it in the end

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
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/Community analyses")

ratsite <- read.table("rats_sites1.csv", sep = ",", head = TRUE)
ratlocs <- ratsite %>% dplyr::select(Longitude, Latitude)
ratlocs <- ratlocs %>% rename(long = Longitude)
ratlocs <- ratlocs %>% rename(lat = Latitude)
rat.xy <- ratlocs %>% as.matrix

# drop 18 sites for functional beta diversity
ratlocs.func <- ratsite %>% filter(Location %ni% c("CMR_Kupe2", "CMR_Kupe3", "ETH_Bale3", "ETH_Bale4", "ETH_Simiens1", "KEN_Turkana2", "LBR_Nimba1", "MLI_Mopti", 
                                             "NAM_Gobabeb", "NAM_Gorrasis", "NAM_Tumasberg", "TZA_KilimanjaroA4", "TZA_KilimanjaroA5", "ZAF_Karoo1", "ZAF_Karoo2", 
                                             "ZAF_Karoo3", "ZAF_Karoo4", "ZAF_Karoo5")) %>% 
  dplyr::select(Longitude, Latitude) %>% rename(long = Longitude) %>% rename(lat = Latitude)
rat.xy.func <- as.matrix(ratlocs.func) 

ratsp1 <-
  read.table("rats_spp1.csv",
             sep = ",",
             head = TRUE,
             row.names = 1)
ratsp <- ifelse(ratsp1 > 0, 1, 0)  # convert to presence-absence

# get environmental (BIOCLIM) data
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/R analysis/Community analyses/R Output/bioclim_data")
enviro.rat1 <- read.csv("bioclim.rat.tax.csv")
enviro.rat <- enviro.rat1[,3:23]

# Select only variables with low VIF (<2.2) based on allsites selection (see Bat_MEM_AM.R)
env.rat.final <- enviro.rat %>% dplyr::select(alt, alt_rough, bio_2, bio_3, 
                                              bio_15, bio_18, bio_19)

# bioclim for functional beta diversity sites (remove 18)
env.rat.func <- enviro.rat1[,3:24] %>% filter(Location %ni% c("CMR_Kupe2", "CMR_Kupe3", "ETH_Bale3", "ETH_Bale4", "ETH_Simiens1", "KEN_Turkana2", "LBR_Nimba1", "MLI_Mopti", 
                                                              "NAM_Gobabeb", "NAM_Gorrasis", "NAM_Tumasberg", "TZA_KilimanjaroA4", "TZA_KilimanjaroA5", "ZAF_Karoo1", "ZAF_Karoo2", 
                                                              "ZAF_Karoo3", "ZAF_Karoo4", "ZAF_Karoo5")) %>% 
  dplyr::select(alt, alt_rough, bio_2, bio_3, bio_15, bio_18, bio_19)

#-------------------------------------------------------------------------------------------------------------------------------------------
# call on some shapefiles for plotting maps

setwd("C:/Dropbox/Ara/GIS_Q/Biogeographic_regions_Africa-Linder")
biogeog1 <- readOGR("Bioregions_Linder_clipped_simplified.shp")
#proj <- "+proj=laea +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
proj1 <- '+proj=longlat + datum=WGS84'
biogeog <- spTransform(biogeog1, CRS(proj1))
proj4string(biogeog) # check crs

#africa1 <- rnaturalearth::ne_countries(continent = "africa", returnclass = "sf")

#-------------------------------------------------------------------------------------------------------------------------------------------
# A few graph-based connectivity schemes (before removing any link):
# ******************************************************************
rat.nbtri <- tri2nb(rat.xy) # Delaunay triangulation
rat.nbgab <- graph2nb(gabrielneigh(rat.xy), sym = TRUE) # Gabriel graph
rat.nbrel <- graph2nb(relativeneigh(rat.xy), sym = TRUE) # Relative neighbourhood graph
rat.nbmst <- mst.nb(dist(rat.xy)) # minimum spanning tree

# Visualisation of the connections:
# ********************************
# The spdep package is designed to create object plotted with base R.
# Visualisation of the connection schemes using ggplot2.
# We first need to transform the nb object into listw objects:
rat.nbtri_listw <- nb2listw(rat.nbtri) 
rat.nbgab_listw <- nb2listw(rat.nbgab)
rat.nbrel_listw <- nb2listw(rat.nbrel)
rat.nbmst_listw <- nb2listw(rat.nbmst)

rat.DA_tri <- nb2ggplot(rat.nbtri_listw, rat.xy)
rat.tri_g <- rat.xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = rat.DA_tri,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "darkred"
  ) +
  labs(title = "Delaunay triangulation") +
  theme_bw()

rat.DA_gab <- nb2ggplot(rat.nbgab_listw, rat.xy)
rat.gab_g <- rat.xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = rat.DA_gab,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "darkred"
  ) +
  labs(title = "Gabriel") +
  theme_bw()

rat.DA_rel <- nb2ggplot(rat.nbrel_listw, rat.xy)
rat.rel_g <- rat.xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = rat.DA_rel,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "chocolate"
  ) +
  labs(title = "Relative neighbourhood") +
  theme_bw()

rat.DA_mst <- nb2ggplot(rat.nbmst_listw, rat.xy)
rat.mst_g <- rat.xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = rat.DA_mst,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "goldenrod"
  ) +
  labs(title = "Min. span tree") +
  theme_bw()
grid.arrange(rat.tri_g, rat.gab_g, rat.rel_g, rat.mst_g, ncol = 2, nrow = 2)

# These commented out lines work in basic R but I didn't change the networks
# nbgab_edited <- edit.nb(nbgab, xy)
# nbrel_edited <- edit.nb(nbrel, xy)
# save(nbgab_edited, nbrel_edited, file = "edited_nb_objects.RData")

# Load the edited nb objects:
# ---------------------------
# load("edited_nb_objects_rats.RData")
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
# nbgab_edited[[2]] # site 2 is connected to sites 18, 88, 101

# Binary forms (no weighths added to the connections):
rat.nbtri_edited_b <- nb2listw(rat.nbtri) # changed names here
rat.nbgab_edited_b <- nb2listw(rat.nbgab) # changed names here
rat.nbrel_edited_b <- nb2listw(rat.nbrel) # changed names here
rat.nbmst_edited_b <- nb2listw(rat.nbmst) # changed names here

#-------------------------------------------------------------------------------------------------------------------------------------------
# Construction of the list of candidate spatial weighting matrices:
# *****************************************************************
# Three candidates:
rat.candidates <- list(rat.tri_b = rat.nbtri_edited_b,
                       rat.gab_b = rat.nbgab_edited_b,
                       rat.rel_b = rat.nbrel_edited_b,
                       rat.nbm_be = rat.nbmst_edited_b)

rat.select <-
  listw.select(
    ratsp,
    rat.candidates,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )

# now save this object 'rat.select' so that we don't have to rerun this!
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/00_Manuscript/00_Journal of Biogeography/Revision/Adam's analyses MEMs")
save(rat.select, file = "rat.select.RData")
load("rat.select.RData")

# Optimised selected SWM:
rat.select$best.id

# Summarises the candidate SWMs:
rat.select$candidates

# Summarises MEM var. selected and Moran's I and p-val in residuals after inclusion of the MEM var.
rat.select$best$summary

# MEM variables to use in models and simulations:
# ***********************************************
rat.MEM.sel <- rat.select$best$MEM.select

# detrend latitude and longitude first
rat.mods <- rda(ratsp ~ long + lat, data = as.data.frame(rat.xy), trace = FALSE)
rat.res <- residuals(rat.mods)

# now run this thru the same analysis from detrended data (above)
rat.select.res <-
  listw.select(
    rat.res,
    rat.candidates,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )
rat.MEM.sel.res <- rat.select.res$best$MEM.select
rat.select.res$candidates
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/00_Manuscript/00_Journal of Biogeography/Revision/Adam's analyses MEMs")
save(rat.select.res, file = "rat.select.res.RData")
load("rat.select.res.RData")

# now run this thru the same analysis but with dudi.pco of beta diversity (rat.beta.pco = taxonomic Bsim)
# beta diversity pco calculated in 'aBetadiversity MEMs_Final.R'
rat.select.beta.tax <-
  listw.select(
    rat.beta.pco$li,
    rat.candidates,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )
rat.MEM.sel.beta.tax <- rat.select.beta.tax$best$MEM.select
rat.select.beta.tax$candidates
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/00_Manuscript/00_Journal of Biogeography/Revision/Adam's analyses MEMs")
save(rat.MEM.sel.beta.tax, file = "rat.MEM.sel.beta.tax.RData")
load("rat.MEM.sel.beta.tax.RData")

# now run this thru the same analysis but with dudi.pco of beta diversity (bats.func.beta.pco = functional Bsim)
# beta diversity pco calculated in 'aBetadiversity MEMs_Final.R'
# but first need a new network because functional diversity dropped two sites
rat.nbtri.func <- nb2listw(tri2nb(rat.xy.func)) # Delaunay triangulation
rat.nbgab.func <- nb2listw(graph2nb(gabrielneigh(rat.xy.func), sym = TRUE)) # Gabriel
rat.nbrel.func <- nb2listw(graph2nb(relativeneigh(rat.xy.func), sym = TRUE)) # Relative neighbourhood
rat.nbmst.func <- nb2listw(mst.nb(dist(rat.xy.func))) # minimum spanning tree

# Four candidates
rat.candidates.func <- list(rat.tri_b = rat.nbtri.func,
                            rat.gab_b = rat.nbgab.func,
                            rat.rel_b = rat.nbrel.func,
                            rat.nbm_b = rat.nbmst.func)

rat.select.beta.func <-
  listw.select(
    rat.beta.func.pco$li,
    rat.candidates.func,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )
rat.MEM.sel.beta.func <- rat.select.beta.func$best$MEM.select
rat.select.beta.func$candidates
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/00_Manuscript/00_Journal of Biogeography/Revision/Adam's analyses MEMs")
save(rat.MEM.sel.beta.func, file = "rat.MEM.sel.beta.func.RData")
load("rat.MEM.sel.beta.func.RData")

# now run this thru the same analysis but with dudi.pco of beta diversity (rat.beta.phyl.pco = phylogenetic Bsim)
# beta diversity pco calculated in 'aBetadiversity MEMs_Final.R'
rat.select.beta.phyl <-
  listw.select(
    rat.beta.phyl.pco$li,
    rat.candidates,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )
rat.MEM.sel.beta.phyl <- rat.select.beta.phyl$best$MEM.select
rat.select.beta.phyl$candidates
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/00_Manuscript/00_Journal of Biogeography/Revision/Adam's analyses MEMs")
save(rat.MEM.sel.beta.phyl, file = "rat.MEM.sel.beta.phyl.RData")
load("rat.MEM.sel.beta.phyl.RData")

#-------------------------------------------------------------------------------------------------------------------------------------------
# Visualisation of the selected MEM variables:
# ********************************************

# s.value() not working, but I've created the MEMs in ggplot (see below)
s.value(
  rat.xy,
  rat.MEM.sel[, c(1:ncol(rat.MEM.sel))],
  ppoint.cex = 0.6,
  symbol = "circle",
  xlim = c(min(rat.xy[, 1]) - 0.3, max(rat.xy[, 1]) + 0.3),
  ylim = c(min(rat.xy[, 2]) - 0.3, max(rat.xy[, 2]) + 0.3)
)

# view in ggplot (taken from 'Moran Eigenvector Map spatial structure.R')
rat.mem.long1 <- rat.MEM.sel %>% mutate('Site' = as.factor(seq(1, nrow(rat.MEM.sel)))) %>% as.data.frame()
rat.mem.long2 <- rat.mem.long1 %>% gather(MEM, Value, MEM1:MEM59)
rat.mem.long2$Longitude <- ratsite$Longitude
rat.mem.long2$Latitude <- ratsite$Latitude

rat.mem.long <- rat.mem.long2 %>% group_by(Longitude, Latitude, MEM) %>% mutate('Mean' = mean(Value)) %>% 
  filter(MEM %in% c("MEM1","MEM6","MEM8","MEM5","MEM7","MEM3","MEM21","MEM4","MEM9",
                    "MEM19","MEM17","MEM16","MEM14","MEM2","MEM15","MEM20","MEM27",
                    "MEM12","MEM11","MEM10","MEM55","MEM13","MEM33","MEM25","MEM18",
                    "MEM32","MEM24","MEM30","MEM29","MEM23","MEM28","MEM26","MEM49",
                    "MEM63","MEM53","MEM42","MEM36","MEM52","MEM37","MEM50","MEM44",
                    "MEM46","MEM60","MEM39","MEM43","MEM22","MEM34","MEM31","MEM40",
                    "MEM56","MEM54","MEM51","MEM57","MEM59"))

rat.mem.long$MEM = factor(rat.mem.long$MEM, levels=c("MEM1","MEM6","MEM8","MEM5","MEM7","MEM3","MEM21","MEM4","MEM9",
                                                     "MEM19","MEM17","MEM16","MEM14","MEM2","MEM15","MEM20","MEM27",
                                                     "MEM12","MEM11","MEM10","MEM55","MEM13","MEM33","MEM25","MEM18",
                                                     "MEM32","MEM24","MEM30","MEM29","MEM23","MEM28","MEM26","MEM49",
                                                     "MEM63","MEM53","MEM42","MEM36","MEM52","MEM37","MEM50","MEM44",
                                                     "MEM46","MEM60","MEM39","MEM43","MEM22","MEM34","MEM31","MEM40",
                                                     "MEM56","MEM54","MEM51","MEM57","MEM59"))

rat.mem.fort <- fortify(rat.mem.long)
biogeog.fort <- merge(fortify(biogeog), as.data.frame(biogeog), by.x="id", by.y=0) # this is needed to add "Region"
#biogeog.fort <- merge(fortify(africa1), as.data.frame(africa1), by.x="sovereignt", by.y=0) # this is needed to add "Region"

rat.mem.graph <- ggplot() +
  # geom_polygon(data=biogeog.fort, aes(x=long, y = lat, group=group), colour='black', 
  #              fill = c('white')) + # I'm failing to give each region a different colour - 'fill = piece' not working
  #coord_fixed(xlim = c(-72,-45), ylim=c(35,70), ratio=1.2)+
  facet_wrap(~MEM, ncol=11, nrow=5) +
  #facet_grid(~MEM) +
  #facet_wrap_paginate(~MEM, page=2) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.title = element_blank())+
  geom_point(data=rat.mem.fort, aes(x = Longitude, y = Latitude, size=Value, fill=Value), shape=22) +
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
rda <- rda(ratsp ~ . , data = MEM.sel)
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
# rat.vp <- varpart(rat.rda$CCA$wa, enviro.rat, rat.MEM.sel) # can use < transfo = "hel" > 
# # but gives identical result, which isn't surprising (Hellinger transformation is for abundance data)
# plot(rat.vp, bg = c(3, 5), Xnames = c("environment", "spatial"))

# repeat but using species-sites data (ie 'batsp') - this should be correct now
rat.vp1 <- varpart(ratsp, env.rat.final, rat.MEM.sel) 
plot(rat.vp1, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))

# now repeat with detrended lat and long (see above) - should lat-long be here? Seems useless
rat.vp1.res <- varpart(rat.res, env.rat.final, rat.MEM.sel.res)
plot(rat.vp1.res, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))

# now repeat with pco of beta diversity (tax =taxonomic, func = functional)
rat.vp1.beta.tax <- varpart(rat.beta.tax.pco$li, env.rat.final, rat.MEM.sel.beta.tax)
rat.vp1.beta.func <- varpart(rat.beta.func.pco$li, env.rat.func, rat.MEM.sel.beta.func)
rat.vp1.beta.phyl <- varpart(rat.beta.phyl.pco$li, env.rat.final, rat.MEM.sel.beta.phyl)

# show rats beta diversity (taxonomic function phylogenetic) side by side
par(mar = c(1.5, 3.5, 1.5, 3.5)) # c(bottom, left, top, right)
par(mfrow = c(3, 1))
plot(rat.vp1.beta.tax, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))
plot(rat.vp1.beta.func, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))
plot(rat.vp1.beta.phyl, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))
par(mfrow = c(1, 1))
par(mar = c(4, 4, 1.5, 0.1)) 


#-------------------------------------------------------------------------------------------------------------------------------------------
# now divide into broad, meso, and fine scale and analyse relative contributions
# of environment and space separately at each scale

# Now repeat but this time use the species sites data (original analysis deleted from above)
# broad scale = MEM: 1, 6, 8, 5 (four)
# meso scale = MEM: 19, 16, 13, 32 (four)
# low scale = MEM: 34, 56, 44, 57 (four)

rat.vp.broad <- varpart(ratsp, env.rat.final, rat.MEM.sel[, c(1,2,3,4)])  
rat.vp.meso <- varpart(ratsp, env.rat.final, rat.MEM.sel[, c(10,12,22,26)])  
rat.vp.fine <- varpart(ratsp, env.rat.final, rat.MEM.sel[, c(47,50,41,53)]) 

par(mar = c(1, 1.5, 1.5, 1)) 
par(mfrow = c(3, 1))
plot(rat.vp.broad, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))
plot(rat.vp.meso, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))
plot(rat.vp.fine, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))
par(mfrow = c(1, 1))
par(mar = c(4, 4, 1.5, 0.1)) 
