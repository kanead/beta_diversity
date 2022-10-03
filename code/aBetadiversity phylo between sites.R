# Phylogenetic betadiversity, code written by Kevin Healy, adapted by Ara Monadjem (15 May 2021)

#------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls()) # removes everything currently held in the R memory
graphics.off() # closes all open graphics windows
# to clear the screen click "control L"
# set.seed(123)

#------------------------------------------------------------------------------------------------------------------------------------------------------

library(ape)
library(caper)
library(phytools)
library(betapart)
library(reshape2)
library(corrplot)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)

#------------------------------------------------------------------------------------------------------------------------------------------------------
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/R analysis/Community analyses/Phylo beta diversity - Kevin code")

#This file contains some of the functions
source("phylostuff.R")

#Read in the phylogeny
full_phylo <- read.tree("cleaned_tree.tre")
species_data <- read.csv("cleaned_species_list.csv", header = T, sep = ",")

#read in presence at site data
bats_presences <- read.csv("Bats_spp1.csv", row.names = 1)
bats_sites <- read.csv("Bats_sites1.csv")

rats_presences <- read.csv("Rats_spp1.csv", row.names = 1)
rats_sites <- read.csv("Rats_sites1.csv")

shrews_presences <- read.csv("Shrews_spp1.csv", row.names = 1)
shrews_sites <- read.csv("Shrews_sites1.csv")

#------------------------------------------------------------------------------------------------------------------------------------------------------
#Match up the inital species list without synonom changes and get the list of unmatched species
phylo_match <- comparative.data(phy = full_phylo,
                                data = species_data,
                                names.col =  "species_data_updated.1")

#Unmatched species
phylo_match$dropped$unmatched.rows

##Synonyms matching
#Lets create a new column called new_name_species to keep track of the change in species names.
#we will create a new species list and overwrite the phylo names below.
new_name_species <- species_data$species_data_updated.1
species_data <- data.frame(species_data, new_name_species)

###Rematch the tree after the synonom changes
species_data_updated <- unique(species_data$new_name_species)

dummy_data <- data.frame(species_data_updated,species_data_updated)

phylo_match_update <- comparative.data(phy = full_phylo,
                                       data = dummy_data,
                                       names.col = "species_data_updated")

phylo_match_update$dropped$unmatched.rows

## Clean species list
cleaned_species_list <- phylo_match_update$data

##clean tree
cleaned_tree <- phylo_match_update$phy

# save this phylogeny for another manuscript
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Hmsc analysis rodents shrews Africa/Data")
ape::write.tree(cleaned_tree, file='small.mammal.phylo.txt') 

## Match the list of species with the cleaned one in the trees
cleaned_bats_presences <- bats_presences[,colnames(bats_presences) %in% cleaned_tree$tip.label]
cleaned_rats_presences <- rats_presences[,colnames(rats_presences) %in% cleaned_tree$tip.label]
cleaned_shrews_presences <- shrews_presences[,colnames(shrews_presences) %in% cleaned_tree$tip.label]

## Combine them into a list for handiness
site_presences <- list(bats = cleaned_bats_presences,
                       rats = cleaned_rats_presences,
                       shrews = cleaned_shrews_presences)

#------------------------------------------------------------------------------------------------------------------------------------------------------
#Now run the actual phylos beta comparisons

# first turn the site data into presence-absence
site_presences$bats[site_presences$bats > "0"] <- c(1)
site_presences$rats[site_presences$rats > "0"] <- c(1)
site_presences$shrews[site_presences$shrews > "0"] <- c(1)

# now run phylo.beta.pair() for phylo beta diversity
bat_phylo_tot_site <- phylo.beta.pair(site_presences$bats, phylo_match_update$phy)
rat_phylo_tot_site <- phylo.beta.pair(site_presences$rats, phylo_match_update$phy)
shrew_phylo_tot_site <- phylo.beta.pair(site_presences$shrews, phylo_match_update$phy)


# to view the individual components of beta diversity for each group
bat.sor1 <- as.matrix(bat_phylo_tot_site$phylo.beta.sor)  
bat.sim1 <- as.matrix(bat_phylo_tot_site$phylo.beta.sim)
bat.nes1 <- as.matrix(bat_phylo_tot_site$phylo.beta.sne)

rat.sor1 <- as.matrix(rat_phylo_tot_site$phylo.beta.sor)  
rat.sim1 <- as.matrix(rat_phylo_tot_site$phylo.beta.sim)
rat.nes1 <- as.matrix(rat_phylo_tot_site$phylo.beta.sne)

shrew.sor1 <- as.matrix(shrew_phylo_tot_site$phylo.beta.sor)  
shrew.sim1 <- as.matrix(shrew_phylo_tot_site$phylo.beta.sim)
shrew.nes1 <- as.matrix(shrew_phylo_tot_site$phylo.beta.sne)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Prepare output as tidy data to write to csv

# bats Bsor (total phylo beta diversity)
bat.psorA <- melt(as.matrix(bat_phylo_tot_site$phylo.beta.sor), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to bat.sorA (above)
batsite.reduced1 <- bats_sites %>%
  dplyr::select(Location, row.Biogeogr = Biogeogr) # using select to also rename biogeogr column
batsite.reduced2 <- bats_sites %>%
  dplyr::select(Location, col.Biogeogr = Biogeogr) # using select to also rename biogeogr column

# now join the biogeographic regions to the phylo betadiversity data.frame bat.fsorA
bat.psorA1 <- bat.psorA %>%
  left_join(batsite.reduced1, by = c('row' = 'Location')) %>%  # batsite.reduced1 and 2 created way above
  left_join(batsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

bat.psor.output <- bat.psorA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
bat.psor.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#bats Bsim (turnover phylo beta diversity)
bat.psimA <- melt(as.matrix(bat_phylo_tot_site$phylo.beta.sim), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to bat.fsimA (done above in bat.sor)
# now join the biogeographic regions to the betadiversity data.frame bat.psimA
bat.psimA1 <- bat.psimA %>%
  left_join(batsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(batsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

bat.psim.output <- bat.psimA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
bat.psim.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#bats Bnes (nestedness phylo beta diversity)
bat.pnesA <- melt(as.matrix(bat_phylo_tot_site$phylo.beta.sne), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to bat.sorA (above)
# now join the biogeographic regions to the betadiversity data.frame bat.nesA
bat.pnesA1 <- bat.pnesA %>%
  left_join(batsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(batsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

bat.pnes.output <- bat.pnesA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
bat.pnes.output %>% print(n = Inf)  

# now create a data.frame with all three betadiversity values in a single data.frame and save to csv for further analysis and graphing
bat.pbetadiversity.sor <- bat.psorA1 %>%
  dplyr::select(row, col, regions, Bsor = value) %>%
  tibble::rownames_to_column("ID")  # need the ID column to allow proper join (since it is a unique value)
bat.pbetadiversity.sim <- bat.psimA1 %>%
  dplyr::select(row, col, regions, Bsim = value) %>%
  tibble::rownames_to_column("ID")
bat.pbetadiversity.sne <- bat.pnesA1 %>%
  dplyr::select(row, col, regions, Bnes = value) %>%
  tibble::rownames_to_column("ID")

bat.pbetadiversity <- bat.pbetadiversity.sor %>%
  inner_join(bat.pbetadiversity.sim, by = "ID") %>%
  inner_join(bat.pbetadiversity.sne, by = "ID") %>%
  dplyr::select(regions, row, col, Bsor, Bsim, Bnes) #%>%
#filter(row != col) # this removes the self-self comparisons (e.g. AGO_Okavango - AGO_Okavango) which will obviously have betadiversity = 0

#------------------------------------------------------------------------------------------------------------------------------------------------------
# rats Bsor (total phylo beta diversity)
rat.psorA <- melt(as.matrix(rat_phylo_tot_site$phylo.beta.sor), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to rat.sorA (above)
ratsite.reduced1 <- rats_sites %>%
  dplyr::select(Location, row.Biogeogr = Biogeogr) # using select to also rename biogeogr column
ratsite.reduced2 <- rats_sites %>%
  dplyr::select(Location, col.Biogeogr = Biogeogr) # using select to also rename biogeogr column

# now join the biogeographic regions to the phylo betadiversity data.frame rat.fsorA
rat.psorA1 <- rat.psorA %>%
  left_join(ratsite.reduced1, by = c('row' = 'Location')) %>%  # ratsite.reduced1 and 2 created way above
  left_join(ratsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

# some comparisons are showing up as NaN; see this by looking a x (below)
x = filter(rat.psorA1, regions == "Congolian-Congolian") 
view(x)

# Hence need to remove NaN from rat.psorA1.
rat.psorA2 <- na.omit(rat.psorA1)

rat.psor.output <- rat.psorA2 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
rat.psor.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#rats Bsim (turnover phylo beta diversity)
rat.psimA <- melt(as.matrix(rat_phylo_tot_site$phylo.beta.sim), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to rat.fsimA (done above in rat.sor)
# now join the biogeographic regions to the betadiversity data.frame rat.psimA
rat.psimA1 <- rat.psimA %>%
  left_join(ratsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(ratsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

# remove NaN from rat.psimA1.
rat.psimA2 <- na.omit(rat.psimA1)

rat.psim.output <- rat.psimA2 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
rat.psim.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#rats Bnes (nestedness phylo beta diversity)
rat.pnesA <- melt(as.matrix(rat_phylo_tot_site$phylo.beta.sne), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to rat.sorA (above)
# now join the biogeographic regions to the betadiversity data.frame rat.nesA
rat.pnesA1 <- rat.pnesA %>%
  left_join(ratsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(ratsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

# remove NaN from rat.pnesA1.
rat.pnesA2 <- na.omit(rat.pnesA1)

rat.pnes.output <- rat.pnesA2 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
rat.pnes.output %>% print(n = Inf)  

# now create a data.frame with all three betadiversity values in a single data.frame and save to csv for further analysis and graphing
rat.pbetadiversity.sor <- na.omit(rat.psorA1) %>%
  dplyr::select(row, col, regions, Bsor = value) %>%
  tibble::rownames_to_column("ID")  # need the ID column to allow proper join (since it is a unique value)
rat.pbetadiversity.sim <- na.omit(rat.psimA1) %>%
  dplyr::select(row, col, regions, Bsim = value) %>%
  tibble::rownames_to_column("ID")
rat.pbetadiversity.sne <- na.omit(rat.pnesA1) %>%
  dplyr::select(row, col, regions, Bnes = value) %>%
  tibble::rownames_to_column("ID")

rat.pbetadiversity <- rat.pbetadiversity.sor %>%
  inner_join(rat.pbetadiversity.sim, by = "ID") %>%
  inner_join(rat.pbetadiversity.sne, by = "ID") %>%
  dplyr::select(regions, row, col, Bsor, Bsim, Bnes) #%>%
#filter(row != col) # this removes the self-self comparisons (e.g. AGO_Okavango - AGO_Okavango) which will obviously have betadiversity = 0

#------------------------------------------------------------------------------------------------------------------------------------------------------
# shrews Bsor (total phylo beta diversity)
shrew.psorA <- melt(as.matrix(shrew_phylo_tot_site$phylo.beta.sor), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to shrew.sorA (above)
shrewsite.reduced1 <- shrews_sites %>%
  dplyr::select(Location, row.Biogeogr = Biogeogr) # using select to also rename biogeogr column
shrewsite.reduced2 <- shrews_sites %>%
  dplyr::select(Location, col.Biogeogr = Biogeogr) # using select to also rename biogeogr column

# now join the biogeographic regions to the phylo betadiversity data.frame shrew.fsorA
shrew.psorA1 <- shrew.psorA %>%
  left_join(shrewsite.reduced1, by = c('row' = 'Location')) %>%  # shrewsite.reduced1 and 2 created way above
  left_join(shrewsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

# remove NaN from shrew.psorA1.
shrew.psorA2 <- na.omit(shrew.psorA1)

shrew.psor.output <- shrew.psorA2 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
shrew.psor.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#shrews Bsim (turnover phylo beta diversity)
shrew.psimA <- melt(as.matrix(shrew_phylo_tot_site$phylo.beta.sim), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to shrew.fsimA (done above in shrew.sor)
# now join the biogeographic regions to the betadiversity data.frame shrew.psimA
shrew.psimA1 <- shrew.psimA %>%
  left_join(shrewsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(shrewsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

# remove NaN from shrew.psimA1.
shrew.psimA2 <- na.omit(shrew.psimA1)

shrew.psim.output <- shrew.psimA2 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
shrew.psim.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#shrews Bnes (nestedness phylo beta diversity)
shrew.pnesA <- melt(as.matrix(shrew_phylo_tot_site$phylo.beta.sne), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to shrew.sorA (above)
# now join the biogeographic regions to the betadiversity data.frame shrew.nesA
shrew.pnesA1 <- shrew.pnesA %>%
  left_join(shrewsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(shrewsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

# remove NaN from shrew.pnesA1.
shrew.pnesA2 <- na.omit(shrew.pnesA1)

shrew.pnes.output <- shrew.pnesA2 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
shrew.pnes.output %>% print(n = Inf)  

# now create a data.frame with all three betadiversity values in a single data.frame and save to csv for further analysis and graphing
shrew.pbetadiversity.sor <- na.omit(shrew.psorA1) %>%
  dplyr::select(row, col, regions, Bsor = value) %>%
  tibble::rownames_to_column("ID")  # need the ID column to allow proper join (since it is a unique value)
shrew.pbetadiversity.sim <- na.omit(shrew.psimA1) %>%
  dplyr::select(row, col, regions, Bsim = value) %>%
  tibble::rownames_to_column("ID")
shrew.pbetadiversity.sne <- na.omit(shrew.pnesA1) %>%
  dplyr::select(row, col, regions, Bnes = value) %>%
  tibble::rownames_to_column("ID")

shrew.pbetadiversity <- shrew.pbetadiversity.sor %>%
  inner_join(shrew.pbetadiversity.sim, by = "ID") %>%
  inner_join(shrew.pbetadiversity.sne, by = "ID") %>%
  dplyr::select(regions, row, col, Bsor, Bsim, Bnes) #%>%
#filter(row != col) # this removes the self-self comparisons (e.g. AGO_Okavango - AGO_Okavango) which will obviously have betadiversity = 0

#------------------------------------------------------------------------------------------------------------------------------------------------------
# save all phylo beta diversity output as csv files

# this combines all components (Bsor, Bsim, Bnes) into a single csv file, separately for each group
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/R analysis/Community analyses/R Output/Betadiversity all_FINAL")

write.csv(bat.pbetadiversity, "bat.phylo.betadiversity.csv")
write.csv(rat.pbetadiversity, "rat.phylo.betadiversity.csv")
write.csv(shrew.pbetadiversity, "shrew.phylo.betadiversity.csv")

setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/R analysis/Community analyses/R Output/Phylo betadiversity - between sites across all regions")

write.csv(bat.psor.output, "bat.phylo_sor.between all sites.csv")
write.csv(bat.psim.output, "bat.phylo_sim.between all sites.csv")
write.csv(bat.pnes.output, "bat.phylo_nes.between all sites.csv")

write.csv(rat.psor.output, "rat.phylo_sor.between all sites.csv")
write.csv(rat.psim.output, "rat.phylo_sim.between all sites.csv")
write.csv(rat.pnes.output, "rat.phylo_nes.between all sites.csv")

write.csv(shrew.psor.output, "shrew.phylo_sor.between all sites.csv")
write.csv(shrew.psim.output, "shrew.phylo_sim.between all sites.csv")
write.csv(shrew.pnes.output, "shrew.phylo_nes.between all sites.csv")

