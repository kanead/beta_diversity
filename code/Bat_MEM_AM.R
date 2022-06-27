# Bat MEM analysis
# code comes from the following paper:
# Benone et al. (2020) How modified landscapes filter rare species and modulate
# the regional pool of ecological traits?

#-------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls()) # removes everything currently held in the R memory
graphics.off() # closes all open graphics windows
# to clear the screen click "control L"
# set.seed(123)

#-------------------------------------------------------------------------------------------------------------------------------------------
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
library(sp)
library(rgdal)
library(sf)  # to project lat/long for variogram
#library(gstat) # to plot variogram

# data(package="vegan")
# data(package="ade4")
# data(package="gstat")


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

batsite <- read.table("Bats_sites1.csv", sep = ",", head = TRUE)
batlocs <- batsite %>% dplyr::select(Longitude, Latitude)
batlocs <- batlocs %>% rename(long = Longitude)
batlocs <- batlocs %>% rename(lat = Latitude)
bat.xy <- batlocs %>% as.matrix

# drop c("ZAR_Kalahari", "ETH_Simien") for functional beta diversity
batlocs.func <- batsite %>% filter(Location %ni% c("ZAR_Kalahari", "ETH_Simien")) %>% 
  dplyr::select(Longitude, Latitude) %>% rename(long = Longitude) %>% rename(lat = Latitude)
bat.xy.func <- as.matrix(batlocs.func) 

batsp1 <-
  read.table("Bats_spp1.csv",
             sep = ",",
             head = TRUE,
             row.names = 1)
batsp <- ifelse(batsp1 > 0, 1, 0)  # convert to presence-absence

# get environmental (BIOCLIM) data
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/R analysis/Community analyses/R Output/bioclim_data")
enviro.bat1 <- read.csv("bioclim.bat.tax.csv")
enviro.bat <- enviro.bat1[,3:23]

# now select those variables based on allsites VIF selection (see below)
env.bat.final <- enviro.bat %>% dplyr::select(alt, alt_rough, bio_2, bio_3, 
                                              bio_15, bio_18, bio_19)

# bioclim for functional beta diversity sites (remove two)
env.bat.func <- enviro.bat1[,3:24] %>% filter(Location %ni% c("ZAR_Kalahari", "ETH_Simien")) %>% 
  dplyr::select(alt, alt_rough, bio_2, bio_3, bio_15, bio_18, bio_19)


# Redo the variable selection based on ALL sites for small mammals (bats, rats, shrews)
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/R analysis/Community analyses/R Output/bioclim_data")
enviro.allsites1 <- read.csv("bioclim.allsites.tax.csv")
enviro.allsites <- enviro.allsites1[,3:23]

vif(enviro.allsites)
# v1 <- vifcor(enviro.allsites, th=0.70) # play around with this t=0.7 means it limits autocorrelation to 0.60
# v2 <- vifstep(enviro.allsites, th=4)  # th=4 here means that VIF is limited to 4 or less

# manually remove variable with highest score but do not delete elevation
# remove bio_11
enviro.allsites2 <- enviro.allsites %>% dplyr::select(-bio_11)
vif(enviro.allsites2)
enviro.allsites3 <- enviro.allsites %>% dplyr::select(-bio_11, - bio_1)
vif(enviro.allsites3)
enviro.allsites4 <- enviro.allsites %>% dplyr::select(-bio_11, - bio_1, -bio_5)
vif(enviro.allsites4)
enviro.allsites5 <- enviro.allsites %>% dplyr::select(-bio_11, - bio_1, -bio_5, -bio_7)
vif(enviro.allsites5)
enviro.allsites6 <- enviro.allsites %>% dplyr::select(-bio_11, - bio_1, -bio_5, -bio_7,
                                                      -bio_6)
vif(enviro.allsites6)
enviro.allsites7 <- enviro.allsites %>% dplyr::select(-bio_11, - bio_1, -bio_5, -bio_7,
                                                      -bio_6, -bio_10)
vif(enviro.allsites7)
enviro.allsites8 <- enviro.allsites %>% dplyr::select(-bio_11, - bio_1, -bio_5, -bio_7,
                                                      -bio_6, -bio_10, -bio_12)
vif(enviro.allsites8)
enviro.allsites9 <- enviro.allsites %>% dplyr::select(-bio_11, - bio_1, -bio_5, -bio_7,
                                                      -bio_6, -bio_10, -bio_12, -bio_17)
vif(enviro.allsites9)
enviro.allsites10 <- enviro.allsites %>% dplyr::select(-bio_11, - bio_1, -bio_5, -bio_7,
                                                      -bio_6, -bio_10, -bio_12, -bio_17,
                                                      -bio_16)
vif(enviro.allsites10)
enviro.allsites11 <- enviro.allsites %>% dplyr::select(-bio_11, - bio_1, -bio_5, -bio_7,
                                                       -bio_6, -bio_10, -bio_12, -bio_17,
                                                       -bio_16, -bio_9)
vif(enviro.allsites11)
enviro.allsites12 <- enviro.allsites %>% dplyr::select(-bio_11, - bio_1, -bio_5, -bio_7,
                                                       -bio_6, -bio_10, -bio_12, -bio_17,
                                                       -bio_16, -bio_9, -bio_4)
vif(enviro.allsites12)
enviro.allsites13 <- enviro.allsites %>% dplyr::select(-bio_11, - bio_1, -bio_5, -bio_7,
                                                       -bio_6, -bio_10, -bio_12, -bio_17,
                                                       -bio_16, -bio_9, -bio_4, -bio_8)
vif(enviro.allsites13)
enviro.allsites14 <- enviro.allsites %>% dplyr::select(-bio_11, - bio_1, -bio_5, -bio_7,
                                                       -bio_6, -bio_10, -bio_12, -bio_17,
                                                       -bio_16, -bio_9, -bio_4, -bio_8,
                                                       -bio_13)
vif(enviro.allsites14)
enviro.allsites15 <- enviro.allsites %>% dplyr::select(-bio_11, - bio_1, -bio_5, -bio_7,
                                                       -bio_6, -bio_10, -bio_12, -bio_17,
                                                       -bio_16, -bio_9, -bio_4, -bio_8,
                                                       -bio_13, -bio_14)
vif(enviro.allsites15)

# the final selection has 7 variables with VIF < 2.2 i.e. no multicollinearity

env.allsites.final <- enviro.allsites %>% dplyr::select(alt, alt_rough, bio_2, bio_3, 
                                           bio_15, bio_18, bio_19)

#-------------------------------------------------------------------------------------------------------------------------------------------
# call on some shapefiles for plotting maps

setwd("C:/Dropbox/Ara/GIS_Q/Biogeographic_regions_Africa-Linder")
biogeog1 <- readOGR("Bioregions_Linder_clipped_simplified.shp")
proj <- "+proj=laea +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
proj1 <- '+proj=longlat + datum=WGS84'
biogeog <- spTransform(biogeog1, CRS(proj1))

bat.loc1 = SpatialPoints(cbind(batsite$Longitude, batsite$Latitude), proj4string=CRS("+proj=longlat"))  # set projection as lat-long
bat.locs <- spTransform(bat.loc1, CRS(proj))
proj4string(biogeog) # check crs

#africa1 <- rnaturalearth::ne_countries(continent = "africa", returnclass = "sf")

#-------------------------------------------------------------------------------------------------------------------------------------------
# A few graph-based connectivity schemes (before removing any link):
# ******************************************************************
bat.nbtri <- tri2nb(bat.xy) # Delaunay triangulation
bat.nbgab <- graph2nb(gabrielneigh(bat.xy), sym = TRUE) # Gabriel
bat.nbrel <- graph2nb(relativeneigh(bat.xy), sym = TRUE) # Relative neighbourhood
bat.nbmst <- mst.nb(dist(bat.xy)) # minimum spanning tree

# Visualisation of the connections:
# ********************************
# The spdep package is designed to create object plotted with base R.
# Visualisation of the connection schemes using ggplot2.
# We first need to transform the nb object into listw objects:
bat.nbtri_listw <- nb2listw(bat.nbtri) 
bat.nbgab_listw <- nb2listw(bat.nbgab)
bat.nbrel_listw <- nb2listw(bat.nbrel)
bat.nbmst_listw <- nb2listw(bat.nbmst)

# create the MEMs based on the network (without species composition data)
bat.mem <- mem(bat.nbmst_listw)
plot(bat.mem)
plot(bat.mem, SpORcoords = bat.locs)
plot(bat.mem[,c(1:18)], SpORcoords = bat.locs)

# note that the Delaunay triangulation is terrible (connections between South & West Africa)
bat.DA_tri <- nb2ggplot(bat.nbtri_listw, bat.xy)
bat.tri_g <- bat.xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = bat.DA_tri,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "darkred"
  ) +
  labs(title = "Delaunay triangulation") +
  theme_bw()

bat.DA_gab <- nb2ggplot(bat.nbgab_listw, bat.xy)
bat.gab_g <- bat.xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = bat.DA_gab,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "darkred"
  ) +
  labs(title = "Gabriel") +
  theme_bw()

bat.DA_rel <- nb2ggplot(bat.nbrel_listw, bat.xy)
bat.rel_g <- bat.xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = bat.DA_rel,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "chocolate"
  ) +
  labs(title = "Relative neighbourhood") +
  theme_bw()

bat.DA_mst <- nb2ggplot(bat.nbmst_listw, bat.xy)
bat.mst_g <- bat.xy %>%
  as.data.frame() %>%
  ggplot(aes(long, lat)) +
  geom_point() +
  coord_fixed() +
  geom_segment(
    data = bat.DA_mst,
    aes(xend = long_to, yend = lat_to),
    size = 0.3,
    alpha = 0.5,
    colour = "goldenrod"
  ) +
  labs(title = "Min. span tree") +
  theme_bw()

grid.arrange(bat.tri_g, bat.gab_g, bat.rel_g, bat.mst_g, ncol = 2, nrow = 2)

# These commented out lines work in basic R but I didn't change the networks
# nbgab_edited <- edit.nb(nbgab, xy)
# nbrel_edited <- edit.nb(nbrel, xy)
# save(nbgab_edited, nbrel_edited, file = "edited_nb_objects.RData")

# Load the edited nb objects:
# ---------------------------
# load("edited_nb_objects_bats.RData")
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
# nbgab_edited[[2]] # site 2 is connected to sites 18, 75, 94, 96

# Binary forms (no weighths added to the connections):
bat.nbtri_edited_b <- nb2listw(bat.nbtri) # changed names here
bat.nbgab_edited_b <- nb2listw(bat.nbgab) # changed names here
bat.nbrel_edited_b <- nb2listw(bat.nbrel) # changed names here
bat.nbmst_edited_b <- nb2listw(bat.nbmst) # changed names here

#-------------------------------------------------------------------------------------------------------------------------------------------
# Construction of the list of candidate spatial weighting matrices:
# *****************************************************************
# Four candidates
bat.candidates <- list(bat.tri_b = bat.nbtri_edited_b,
                       bat.gab_b = bat.nbgab_edited_b,
                       bat.rel_b = bat.nbrel_edited_b,
                       bat.nbm_b = bat.nbmst_edited_b)

# Optimisation of a subset of spatial predictors
# optimise the selection of a subset of spatial predictors from the
# best-suited SWM.
bat.select <-
  listw.select(
    batsp,
    bat.candidates,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )

# now save this object 'bat.select' so that we don't have to rerun this!
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/00_Manuscript/00_Journal of Biogeography/Revision/Adam's analyses MEMs")
save(bat.select, file = "bat.select.RData")
load("bat.select.RData")

# Optimised selected SWM:
bat.select$best.id

# Summarises the candidate SWMs:
bat.select$candidates

# Summarises MEM var. selected and Moran's I and p-val in residuals after inclusion of the MEM var.
bat.select$best$summary

# MEM variables to use in models and simulations:
# ***********************************************
bat.MEM.sel <- bat.select$best$MEM.select

# detrend latitude and longitude first
bat.mods <- rda(batsp ~ long + lat, data = as.data.frame(bat.xy), trace = FALSE)
bat.res <- residuals(bat.mods)

# now run this thru the same analysis but with detrended data (from above)
bat.select.res <-
  listw.select(
    bat.res,
    bat.candidates,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )
bat.MEM.sel.res <- bat.select.res$best$MEM.select
bat.select.res$candidates
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/00_Manuscript/00_Journal of Biogeography/Revision/Adam's analyses MEMs")
save(bat.select.res, file = "bat.select.res.RData")
load("bat.select.res.RData")

# now run this thru the same analysis but with dudi.pco of beta diversity (bat.beta.pco = taxonomic Bsim)
# beta diversity pco calculated in 'aBetadiversity MEMs_Final.R'
bat.select.beta.tax <-
  listw.select(
    bat.beta.tax.pco$li,
    bat.candidates,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )
bat.MEM.sel.beta.tax <- bat.select.beta.tax$best$MEM.select
bat.select.beta.tax$candidates
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/00_Manuscript/00_Journal of Biogeography/Revision/Adam's analyses MEMs")
save(bat.MEM.sel.beta.tax, file = "bat.MEM.sel.beta.tax.RData")
load("bat.MEM.sel.beta.tax.RData")

# now run this thru the same analysis but with dudi.pco of beta diversity (bats.func.beta.pco = functional Bsim)
# beta diversity pco calculated in 'aBetadiversity MEMs_Final.R'
# but first need a new network because functional diversity dropped two sites

bat.nbtri.func <- nb2listw(tri2nb(bat.xy.func)) # Delaunay triangulation
bat.nbgab.func <- nb2listw(graph2nb(gabrielneigh(bat.xy.func), sym = TRUE)) # Gabriel
bat.nbrel.func <- nb2listw(graph2nb(relativeneigh(bat.xy.func), sym = TRUE)) # Relative neighbourhood
bat.nbmst.func <- nb2listw(mst.nb(dist(bat.xy.func))) # minimum spanning tree

# Four candidates
bat.candidates.func <- list(bat.tri_b = bat.nbtri.func,
                       bat.gab_b = bat.nbgab.func,
                       bat.rel_b = bat.nbrel.func,
                       bat.nbm_b = bat.nbmst.func)

bat.select.beta.func <-
  listw.select(
    bat.beta.func.pco$li,
    bat.candidates.func,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )
bat.MEM.sel.beta.func <- bat.select.beta.func$best$MEM.select
bat.select.beta.func$candidates
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/00_Manuscript/00_Journal of Biogeography/Revision/Adam's analyses MEMs")
save(bat.MEM.sel.beta.func, file = "bat.MEM.sel.beta.func.RData")
load("bat.MEM.sel.beta.func.RData")

# now run this thru the same analysis but with dudi.pco of beta diversity (bat.beta.phyl.pco = phylogenetic Bsim)
# beta diversity pco calculated in 'aBetadiversity MEMs_Final.R'
bat.select.beta.phyl <-
  listw.select(
    bat.beta.phyl.pco$li,
    bat.candidates,
    MEM.autocor = "positive",
    method = "FWD",
    p.adjust = TRUE
  )
bat.MEM.sel.beta.phyl <- bat.select.beta.phyl$best$MEM.select
bat.select.beta.phyl$candidates
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/00_Manuscript/00_Journal of Biogeography/Revision/Adam's analyses MEMs")
save(bat.MEM.sel.beta.phyl, file = "bat.MEM.sel.beta.phyl.RData")
load("bat.MEM.sel.beta.phyl.RData")


#-------------------------------------------------------------------------------------------------------------------------------------------
# Visualisation of the selected MEM variables: [s.value() is not working anymore!]
# ********************************************

# s.value() not working, but I've created the MEMs in ggplot (see below)
s.value(
  bat.xy,
  bat.MEM.sel[, c(1:ncol(bat.MEM.sel))],
#  ppoint.cex = 0.6,
#  symbol = "circle",
  xlim = c(min(bat.xy[, 1]) - 0.3, max(bat.xy[, 1]) + 0.3),
  ylim = c(min(bat.xy[, 2]) - 0.3, max(bat.xy[, 2]) + 0.3)
)

# view in ggplot (taken from 'Moran Eigenvector Map spatial structure.R')
bat.mem.long1 <- bat.MEM.sel %>% mutate('Site' = as.factor(seq(1, nrow(bat.MEM.sel)))) %>% as.data.frame()
bat.mem.long2 <- bat.mem.long1 %>% gather(MEM, Value, MEM1:MEM23)
bat.mem.long2$Longitude <- batsite$Longitude
bat.mem.long2$Latitude <- batsite$Latitude

bat.mem.long <- bat.mem.long2 %>% group_by(Longitude, Latitude, MEM) %>% mutate('Mean' = mean(Value)) %>% 
  filter(MEM %in% c("MEM1","MEM2","MEM6","MEM3","MEM8","MEM4","MEM5","MEM10","MEM9","MEM20",
                    "MEM14","MEM13","MEM19","MEM24","MEM7","MEM12","MEM15","MEM16","MEM26",
                    "MEM11","MEM18","MEM21","MEM22","MEM27","MEM23"))

bat.mem.long$MEM = factor(bat.mem.long$MEM, levels=c("MEM1","MEM2","MEM6","MEM3","MEM8","MEM4","MEM5","MEM10","MEM9","MEM20",
                                                     "MEM14","MEM13","MEM19","MEM24","MEM7","MEM12","MEM15","MEM16","MEM26",
                                                     "MEM11","MEM18","MEM21","MEM22","MEM27","MEM23"))

bat.mem.fort <- fortify(bat.mem.long)
biogeog.fort <- merge(fortify(biogeog), as.data.frame(biogeog), by.x="id", by.y=0) # this is needed to add "Region"
#biogeog.fort <- merge(fortify(africa1), as.data.frame(africa1), by.x="sovereignt", by.y=0) # this is needed to add "Region"

bat.mem.graph <- ggplot() +
  # geom_polygon(data=biogeog.fort, aes(x=long, y = lat, group=group), colour='black', 
  #              fill = c('white')) + # I'm failing to give each region a different colour - 'fill = piece' not working
  #coord_fixed(xlim = c(-72,-45), ylim=c(35,70), ratio=1.2)+
  facet_wrap(~MEM, ncol=10, nrow=3) +
  theme(panel.background = element_rect(fill = "white", colour = "black"),
        axis.title = element_blank())+
  geom_point(data=bat.mem.fort, aes(x = Longitude, y = Latitude, size=Value, fill=Value), shape=22) +
  theme_bw()+theme(legend.position = "right",
                   panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(y= 'Latitude') +  
  labs(x= 'Longitude') +
  theme(axis.text.x=element_text(colour="white")) +
  theme(axis.text.y=element_text(colour="white")) +
  scale_fill_gradient2(low='black', mid='white', high="red", midpoint = 0) #+ # Adam, how can I manually set the size of the dots based on the value of 'Value'? I want to exaggerate the size difference

#-------------------------------------------------------------------------------------------------------------------------------------------
# use variogram to estimate the relative distance of each MEM
# first convert lat-long to projected coordinates - Lambert-Zimuthal Equal Area

bat.xy.proj <- as.data.frame(bat.xy) %>% st_as_sf(coords = c('long', 'lat'), crs = 4326) %>%
  st_transform("+proj=laea +x_0=0 +y_0=0 +lon_0=0 +lat_0=0")

# now extract lat and long from geometry
bat.xy.proj1 <- as.data.frame(st_coordinates(bat.xy.proj$geometry)) 

# examples of variograms done manually
bat.MEM1.vario <- variogmultiv(bat.MEM.sel[,1], bat.xy.proj1, dmax = max(dist(bat.xy.proj1)), 
                          nclass = 20)
bat.MEM10.vario <- variogmultiv(bat.MEM.sel[,8], bat.xy.proj1, dmax = max(dist(bat.xy.proj1)), 
                               nclass = 20)
bat.MEM22.vario <- variogmultiv(bat.MEM.sel[,24], bat.xy.proj1, dmax = max(dist(bat.xy.proj1)), 
                                nclass = 20)
bat.MEM27.vario <- variogmultiv(bat.MEM.sel[,23], bat.xy.proj1, dmax = max(dist(bat.xy.proj1)), 
                                nclass = 20)

options(scipen=5)  # set scipen to avoid scientific notation
par(mar = c(1, 1, 1.5, 1.5))
par(mfrow = c(4, 1))
plot(bat.MEM1.vario$d/1000, bat.MEM1.vario$var, pch = 20, 
     xlab = "Distance (km)", ylab = "Semivariance", ty = 'b', main = 'broad scale')
plot( bat.MEM10.vario$d/1000, bat.MEM10.vario$var, ty = 'b', pch = 20, 
      xlab = "Distance (km)", ylab = "Semivariance", main = 'meso scale')
plot( bat.MEM27.vario$d/1000, bat.MEM27.vario$var, ty = 'b', pch = 20, 
      xlab = "Distance (km)", ylab = "Semivariance", main = 'fine scale')
par(mfrow = c(1, 1))
par(mar = c(4, 4, 1.5, 0.1))   

# now run variograms in for loop, based on taxonomic beta diversity (i.e. bat.MEM.sel.beta.tax)
bat.tax.vario <- list()
#values <- c(3,2,1,4,9,5,11,8,13,7,6,16,12,25,39,24,18,19,29,21,10)  # these are the names of the 22 significant MEMs from bat.MEM.sel.beta.tax
values <- c(1:21)  # these are the col numbers of the 21 significant MEMs from bat.MEM.sel.beta.tax

for (i in values) {
  output <- variogmultiv(bat.MEM.sel.beta.tax[,i], bat.xy.proj1, dmax = max(dist(bat.xy.proj1)),
                         nclass = 10)
  bat.tax.vario[[i]] <- output        # Store output in list
}

bat.tax.vario
#bat.tax.vario %>% discard(is.null)

par(mar = c(1, 1, 1.5, 1.5))
par(mfrow = c(5, 5))

plot(bat.tax.vario[[3]][['d']]/1000, bat.tax.vario[[3]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM3", ty = 'b')
plot(bat.tax.vario[[2]][['d']]/1000, bat.tax.vario[[2]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM2", ty = 'b')
plot(bat.tax.vario[[1]][['d']]/1000, bat.tax.vario[[1]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM1", ty = 'b')
plot(bat.tax.vario[[4]][['d']]/1000, bat.tax.vario[[4]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM4", ty = 'b')
plot(bat.tax.vario[[9]][['d']]/1000, bat.tax.vario[[9]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM9", ty = 'b')
plot(bat.tax.vario[[5]][['d']]/1000, bat.tax.vario[[5]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM5", ty = 'b')
plot(bat.tax.vario[[11]][['d']]/1000, bat.tax.vario[[11]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM11", ty = 'b')
plot(bat.tax.vario[[8]][['d']]/1000, bat.tax.vario[[8]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM8", ty = 'b')
plot(bat.tax.vario[[13]][['d']]/1000, bat.tax.vario[[13]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM13", ty = 'b')
plot(bat.tax.vario[[7]][['d']]/1000, bat.tax.vario[[7]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM7", ty = 'b')
plot(bat.tax.vario[[6]][['d']]/1000, bat.tax.vario[[6]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM6", ty = 'b')
plot(bat.tax.vario[[16]][['d']]/1000, bat.tax.vario[[16]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM16", ty = 'b')
plot(bat.tax.vario[[12]][['d']]/1000, bat.tax.vario[[12]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM12", ty = 'b')
plot(bat.tax.vario[[25]][['d']]/1000, bat.tax.vario[[25]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM25", ty = 'b')
plot(bat.tax.vario[[39]][['d']]/1000, bat.tax.vario[[39]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM39", ty = 'b')
plot(bat.tax.vario[[24]][['d']]/1000, bat.tax.vario[[24]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM24", ty = 'b')
plot(bat.tax.vario[[18]][['d']]/1000, bat.tax.vario[[18]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM18", ty = 'b')
plot(bat.tax.vario[[19]][['d']]/1000, bat.tax.vario[[19]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM19", ty = 'b')
plot(bat.tax.vario[[29]][['d']]/1000, bat.tax.vario[[29]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM29", ty = 'b')
plot(bat.tax.vario[[21]][['d']]/1000, bat.tax.vario[[21]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM21", ty = 'b')
plot(bat.tax.vario[[10]][['d']]/1000, bat.tax.vario[[10]][['var']], pch = 20,
     xlab = "Distance (km)", ylab = "MEM10", ty = 'b')

par(mfrow = c(1, 1))
par(mar = c(4, 4, 1.5, 0.1))   


#-------------------------------------------------------------------------------------------------------------------------------------------
# this section of correlograms is not working and can be ignored or deleted
# now run correlograms (not sure when to use correlogram and when variogram!)

# redo the network analysis so as to define a threshold
bat.net <- dnearneigh(bat.xy, 0, 0.1, longlat=TRUE)  # all links between 0 and 650km

# bat.net1 <- dnearneigh(bat.xy.proj, 0, 650,  # all links between 0 and 650km
#                        use_s2=FALSE, sf::st_is_longlat(x) == TRUE, sf::sf_use_s2() == TRUE) 



bat.MEM1.corr <- sp.correlogram(bat.net, bat.MEM.sel$MEM1, order=10, method = "I",
                                style = "W", randomisation=FALSE)
print(bat.MEM1.corr, p.adj.method = "holm")
plot(bat.MEM1.corr)


proj <- "+proj=laea +lon_0=0 +lat_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs"
bat.loc1 <- st_as_sf(x=batlocs,                          
                     coords = c("Longitude", "Latitude"),
                     crs = "+proj=longlat +ellps=WGS84 +no_defs")

bat.loc2 <- st_transform(bat.loc1, crs = proj)

# bat.net1 <- graph2nb(gabrielneigh(bat.xy), sym=TRUE)

# networks are refusing to plot so I can't see why they are not working above
s.label(bat.loc1$geometry, nb=bat.net)

#-------------------------------------------------------------------------------------------------------------------------------------------
# Visualisation of the species composition multiscale spatial patterns:
# *********************************************************************
bat.rda.mem.up <- rda(batsp ~ ., data = bat.MEM.sel)
bat.rda.mem.lw <- rda(batsp ~ 1, data = bat.MEM.sel)
set.seed(1)
bat.mods <- ordistep(bat.rda.mem.lw, scope = formula(bat.rda.mem.up),
                     trace = FALSE)
bat.mods$anova

# apparently the stepwise selection (above) is not appropriate
# need a 2-step selection process and hence ordiR2step() - see below
bat.mods1 <- ordiR2step(bat.rda.mem.lw, bat.rda.mem.up, trace = FALSE)
bat.mods1$anova

# now run rda on environmental variables
bat.rda.env.lw <- rda(batsp ~ 1, data = env.bat.final)
bat.rda.env.up <- rda(batsp ~ ., data = env.bat.final)
set.seed(1)
bat.mods.env <- ordiR2step(bat.rda.env.lw, bat.rda.env.up, trace = FALSE)
bat.mods.env$anova # not all 6 bioclim variables are significant (bioclim14 dropped)

# what does it look like running rda on all bioclim variables?
bat.rda.env.lw1 <- rda(batsp ~ 1, data = enviro.bat)
bat.rda.env.up1 <- rda(batsp ~ ., data = enviro.bat)
set.seed(1)
bat.mods.env1 <- ordiR2step(bat.rda.env.lw1, bat.rda.env.up1, trace = FALSE)
bat.mods.env1$anova

# Interesting to compare bat.mods.env with bat.mods.env1 (bio2 and 3 are dropped and alt gained!)

# can also run rda on environmental data
env.rda <- rda(env.bat.final ~ . , data = MEM.sel) # and can repeat this for the 3 scales

# subset for each of the 3 spatial scales
env.rda <- rda(env.bat.final ~ . , data = MEM.sel[,1:4])

#-------------------------------------------------------------------------------------------------------------------------------------------
# Permutation test (although the correct p-value is the adjusted one obtained with listw.select)
anova.cca(bat.rda)

# Permutation test of each constraint axis (RDA axes):
test_by_axis <- anova.cca(bat.rda, by = "axis")
nb_ax <- length(which(test_by_axis$"Pr(>F)" <= 0.05))
paste("The first", nb_ax, "RDA axes are significant.", sep = " ")

# Adjusted R-squared of the model:
RsquareAdj(bat.rda)$adj.r.squared # must correspond to AdjR2 obtained with select$best$summary
#' eigenvalues of RDAs 
sc <- scores(
  bat.rda,
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
# repeat but using species-sites data (ie 'batsp') - this should be correct now
bat.vp1 <- varpart(batsp, env.bat.final, bat.MEM.sel)
plot(bat.vp1, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))

# now repeat with detrended lat and long (see above) - should lat-long be here? Seems useless
bat.vp1.res <- varpart(bat.res, env.bat.final, bat.MEM.sel.res)
plot(bat.vp1.res, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))

# now repeat with pco of beta diversity (tax =taxonomic, func = functional)
# bat.beta.tax.pco and all other .pco calculated in 'aBetadiversity MEMs_Final.R'
bat.vp1.beta.tax <- varpart(bat.beta.tax.pco$li, env.bat.final, bat.MEM.sel.beta.tax)
bat.vp1.beta.func <- varpart(bat.beta.func.pco$li, env.bat.func, bat.MEM.sel.beta.func)
bat.vp1.beta.phyl <- varpart(bat.beta.phyl.pco$li, env.bat.final, bat.MEM.sel.beta.phyl)

# show bats beta diversity (taxonomic function phylogenetic) side by side
par(mar = c(1.5, 3.5, 1.5, 3.5)) # c(bottom, left, top, right)
par(mfrow = c(3, 1))
plot(bat.vp1.beta.tax, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))
plot(bat.vp1.beta.func, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))
plot(bat.vp1.beta.phyl, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))
par(mfrow = c(1, 1))
par(mar = c(4, 4, 1.5, 0.1)) 

# show bats, rats, shrews beta diversity (taxonomic and function) side by side
par(mar = c(1, 1.5, 1.5, 1)) 
par(mfrow = c(3, 2))
plot(bat.vp1.beta.tax, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))
plot(bat.vp1.beta.func, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))
plot(rat.vp1.beta.tax, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))
plot(rat.vp1.beta.func, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))
plot(shrew.vp1.beta.tax, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))
plot(shrew.vp1.beta.func, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial", "lat-long"))
par(mfrow = c(1, 1))
par(mar = c(4, 4, 1.5, 0.1)) 




# # we get a similar result if we use all bioclim variables (ie with collinearity)
# bat.vp2 <- varpart(batsp, enviro.bat, bat.MEM.sel)
# plot(bat.vp2, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))

# # separate elevation from bioclim variables as a 3rd partition (not useful as alt doesn't explain much)
# alt.bat <- env.bat.final[,1:2]
# bioclim.bat <- env.bat.final[,3:7]
# bat.vp3 <- varpart(batsp, alt.bat, bioclim.bat, bat.MEM.sel)
# #showvarparts(3)
# plot(bat.vp3, bg = c('grey5', 'grey45', 'grey85'), Xnames = c('alt', 'bioclim', 'space'))

# plot bats, rats and shrews side by side
par(mar = c(1, 1, 1.5, 1.5)) 
par(mfrow = c(3, 1))
plot(bat.vp1, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))
plot(rat.vp1, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))
plot(shrew.vp1, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))
par(mfrow = c(1, 1))
par(mar = c(4, 4, 1.5, 0.1)) 

# plot bats, rats and shrews side by side (but this time with detrended lat-long)
par(mar = c(1, 1, 1.5, 1.5)) 
par(mfrow = c(3, 1))
plot(bat.vp1.res, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))
plot(rat.vp1.res, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))
plot(shrew.vp1.res, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))
par(mfrow = c(1, 1))
par(mar = c(4, 4, 1.5, 0.1)) 

#-------------------------------------------------------------------------------------------------------------------------------------------
# now divide into broad, meso, and fine scale and analyse relative contributions
# of environment and space separately at each scale

# I think the above is not correct (deleted). Now repeat but this time use the species sites data
# base scale on MEMs (ie. visualizing 'bat.MEM.sel')

# broad scale = MEM: 1, 10, 11, 12 (four)
# meso scale = MEM: 13, 14, 15, 6 (four)
# low scale = MEM: 16, 18, 20, 24 (four)

bat.vp.broad <- varpart(batsp, env.bat.final, bat.MEM.sel[, c(1,8,16,20)])
bat.vp.meso <- varpart(batsp, env.bat.final, bat.MEM.sel[, c(12,11,17,3)])
bat.vp.fine <- varpart(batsp, env.bat.final, bat.MEM.sel[, c(18,21,10,14)])

par(mar = c(1, 1.5, 1.5, 1)) 
par(mfrow = c(3, 1))
plot(bat.vp.broad, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))
plot(bat.vp.meso, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))
plot(bat.vp.fine, bg = c('grey5', 'grey95'), Xnames = c("environment", "spatial"))
par(mfrow = c(1, 1))
par(mar = c(4, 4, 1.5, 0.1)) 
