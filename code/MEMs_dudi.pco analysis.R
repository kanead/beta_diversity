# Testing taxonomic, functional and phylogenetic beta diversity using Moran's Eigenvector Maps - for the revision of the J Biogeogr manuscript
# Specifically, for running the dudi.pco (principal coordinates analysis) that then feeds into the listw() function in the MEM_bat_AM, MEM_rat_AM and MEM_shrew_AM


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
#library(sp)
library(rgdal)
library(sf)  # to project lat/long for variogram

#-------------------------------------------------------------------------------------------------------------------------------------------
# load data
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/Community analyses")

batsite <- read.table("Bats_sites1.csv", sep = ",", head = TRUE)
batlocs <- batsite %>% dplyr::select(Longitude, Latitude)

batsp1 <- read.csv("Bats_spp1.csv", row.names = 1)
batsp <- ifelse(batsp1 > 0, 1, 0)  # convert to presence-absence

setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/R analysis/Community analyses/R Output/Betadiversity all_FINAL")

betadiversity <- read.csv("betadiversity_all.csv", row.names = 1)

bat.tax <- read.csv("bat.taxonomic.betadiversity.csv")
rat.tax <- read.csv("rat.taxonomic.betadiversity.csv")
shrew.tax <- read.csv("shrew.taxonomic.betadiversity.csv")

bat.func <- read.csv("bat.functional.betadiversity.csv")
rat.func <- read.csv("rat.functional.betadiversity.csv")
shrew.func <- read.csv("shrew.functional.betadiversity.csv")

bat.phylo <- read.csv("bat.phylo.betadiversity.csv")
rat.phylo <- read.csv("rat.phylo.betadiversity.csv")
shrew.phylo <- read.csv("shrew.phylo.betadiversity.csv")



setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/00_Manuscript/00_Journal of Biogeography/Revision1/Adam's analyses MEMs")
load("bat.select.res.RData")
load("rat.select.res.RData")
load("shrew.select.res.RData")

bat.select.res

#-------------------------------------------------------------------------------------------------------------------------------------------
# obtaining the matrix of beta diversity (from betapart) - need to run this first separately ('aBetadiversity between sites_FINAL.R')
# - and also need to run 'aBetadiversity phylo between sites_Based on Kevin_FINAL.R'

bat.beta.tax # taxonomic beta diversity of bats
rat.beta.tax # taxonomic beta diversity of rats
shrew.beta.tax # taxonomic beta diversity of shrews

bat.beta.func # functional beta diversity of bats
rat.beta.func # functional beta diversity of rats
shrew.beta.func # functional beta diversity of shrews

bat.beta.phyl <- bat_phylo_tot_site # phylogenetic beta diversity of bats
rat.beta.phyl <- rat_phylo_tot_site # phylogenetic beta diversity of rats
shrew.beta.phyl <- shrew_phylo_tot_site # phylogenetic beta diversity of shrews

# now run dudi.pco (convert to principal coordinates see Garcia-Giron et al (2019) Science of the Total Environment)
# for taxonomic beta diversity
bat.beta.tax.pco <-dudi.pco(bat.beta.tax$beta.sim, scannf=FALSE, nf=10)
rat.beta.tax.pco <-dudi.pco(rat.beta.tax$beta.sim, scannf=FALSE, nf=10)
shrew.beta.tax.pco <-dudi.pco(shrew.beta.tax$beta.sim, scannf=FALSE, nf=10)

# now run for functional beta diversity
bat.beta.func.pco <- dudi.pco(bat.beta.func$funct.beta.sim, scannf=FALSE, nf=10)
rat.beta.func.pco <- dudi.pco(rat.beta.func$funct.beta.sim, scannf=FALSE, nf=10)
shrew.beta.func.pco <- dudi.pco(shrew.beta.func$funct.beta.sim, scannf=FALSE, nf=10)

# now run for phylogenetic beta diversity
bat.beta.phyl.pco <- dudi.pco(bat.beta.phyl$phylo.beta.sim, scannf=FALSE, nf=10)
rat.beta.phyl.pco <- dudi.pco(rat.beta.phylaa, scannf=FALSE, nf=10)
shrew.beta.phyl.pco <- dudi.pco(shrew.beta.phylaa, scannf=FALSE, nf=10)

# sorting out the error with missing value in d in above analysis (due to NaNs)
# first convert matrix to dataframe using melt
rat.beta.phyl1 <- melt(as.matrix(rat.beta.phyl$phylo.beta.sim), varnames = c("row", "col")) # notice the NaNs!
mean(na.omit(rat.beta.phyl1$value)) # calculate the mean value
rat.beta.phyl1[is.na(rat.beta.phyl1)] <- 0.3437805 # impute mean value for NaNs

rat.beta.phyl0 <- dcast(rat.beta.phyl1, row~col) %>% column_to_rownames('row') # rearrange dataframe to wide format
rat.beta.phylaa <- as.dist(rat.beta.phyl0) # and convert to class 'dist'

# and now for shrews
shrew.beta.phyl1 <- melt(as.matrix(shrew.beta.phyl$phylo.beta.sim), varnames = c("row", "col")) # notice the NaNs!
mean(na.omit(shrew.beta.phyl1$value))
shrew.beta.phyl1[is.na(shrew.beta.phyl1)] <- 0.2468167 # impute mean value for NaNs

shrew.beta.phyl0 <- dcast(shrew.beta.phyl1, row~col) %>% column_to_rownames('row')
shrew.beta.phylaa <- as.dist(shrew.beta.phyl0)
