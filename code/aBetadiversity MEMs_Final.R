# Testing taxonomic, functional and phylogenetic beta diversity using Moran's Eigenvector Maps

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



setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/00_Manuscript/00_Journal of Biogeography/Revision/Adam's analyses MEMs")
load("bat.select.res.RData")
load("rat.select.res.RData")
load("shrew.select.res.RData")

bat.select.res

#-------------------------------------------------------------------------------------------------------------------------------------------
# obtaining the matrix of beta diversity (from betapart) - need to run this first separately ('aBetadiversity between sites_FINAL.R')

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

# trying to figure out the error with missing value in d
rat.beta.phyl1 <- melt(as.matrix(rat.beta.phyl$phylo.beta.sim), varnames = c("row", "col")) # notice the NaNs!
rat.beta.phyl1[is.na(rat.beta.phyl1)] <- 0.5 # impute 0.5 for NaNs

rat.beta.phyl0 <- dcast(rat.beta.phyl1, row~col) %>% column_to_rownames('row')
rat.beta.phylaa <- as.dist(rat.beta.phyl0)

# and now for shrews
shrew.beta.phyl1 <- melt(as.matrix(shrew.beta.phyl$phylo.beta.sim), varnames = c("row", "col")) # notice the NaNs!
shrew.beta.phyl1[is.na(shrew.beta.phyl1)] <- 0.5 # impute 0.5 for NaNs

shrew.beta.phyl0 <- dcast(shrew.beta.phyl1, row~col) %>% column_to_rownames('row')
shrew.beta.phylaa <- as.dist(shrew.beta.phyl0)
