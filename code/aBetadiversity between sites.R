# Beta diversity between groups using beta.pair() i.e. site by site comparisons, Ara Monadjem (12 May 2021)

#------------------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls()) # removes everything currently held in the R memory
graphics.off() # closes all open graphics windows
# to clear the screen click "control L"
# set.seed(123)
#------------------------------------------------------------------------------------------------------------------------------------------------------

library(betapart)
library(ecodist)
library(vegan)
library(tidyverse)
library(sf)
library(RVAideMemoire)  # not working properly (for me)
#library(pairwise.adonis) # not working properly (for me)
library(reshape2)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# input data community data

setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/Community analyses")

batsp <- read.table("Bats_spp1.csv", sep=",", head = TRUE, row.names = 1)
bat_pres_abs <- ifelse(batsp>0,1,0)  # convert to presence-absence
batsite <- read.table("Bats_sites1.csv", sep=",", head=TRUE)
bats.traits <- read.table("Bats_traits1.csv", sep=",", head=TRUE)

ratsp <- read.table("Rats_spp1.csv", sep=",", head = TRUE, row.names = 1)
rat_pres_abs <- ifelse(ratsp>0,1,0)
ratsite <- read.table("Rats_sites1.csv", sep=",", head = TRUE)
rats.traits <- read.table("Rats_traits1.csv", sep=",", head=TRUE)

shrewsp <- read.table("Shrews_spp1.csv", sep=",", head = TRUE, row.names = 1)
shrew_pres_abs <- ifelse(shrewsp>0,1,0)
shrewsite <- read.table("Shrews_sites1.csv", sep=",", head = TRUE)
shrews.traits <- read.table("Shrews_traits1.csv", sep=",", head=TRUE)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Beta diversity calculations using betapart beta.pair(), which compares between all sites (i.e. does not combine sites into a single output value)

bat.beta.tax <- beta.pair(bat_pres_abs, index.family="sorensen")
rat.beta.tax <- beta.pair(rat_pres_abs, index.family="sorensen")
shrew.beta.tax <- beta.pair(shrew_pres_abs, index.family="sorensen")

#------------------------------------------------------------------------------------------------------------------------------------------------------
#bats Bsor (total beta diversity)
bat.sorA <- melt(as.matrix(bat.beta$beta.sor), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to bat.sorA (above)
batsite.reduced1 <- batsite %>%
  dplyr::select(Location, row.Biogeogr = Biogeogr) # using select to also rename biogeogr column
batsite.reduced2 <- batsite %>%
  dplyr::select(Location, col.Biogeogr = Biogeogr) # using select to also rename biogeogr column

# now join the biogeographic regions to the betadiversity data.frame bat.sorA
bat.sorA1 <- bat.sorA %>%
  left_join(batsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(batsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

bat.sor.output <- bat.sorA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
bat.sor.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#bats Bsim (turnover beta diversity)
bat.simA <- melt(as.matrix(bat.beta$beta.sim), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to bat.simA (done above in bat.sor)
# now join the biogeographic regions to the betadiversity data.frame bat.simA
bat.simA1 <- bat.simA %>%
  left_join(batsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(batsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

bat.sim.output <- bat.simA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
bat.sim.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#bats Bnes (nestedness beta diversity)
bat.nesA <- melt(as.matrix(bat.beta$beta.sne), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to bat.sorA (above)
# now join the biogeographic regions to the betadiversity data.frame bat.nesA
bat.nesA1 <- bat.nesA %>%
  left_join(batsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(batsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

bat.nes.output <- bat.nesA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
bat.nes.output %>% print(n = Inf)  

# now create a data.frame with all three betadiversity values in a single data.frame and save to csv for further analysis and graphing
bat.betadiversity.sor <- bat.sorA1 %>%
  dplyr::select(row, col, regions, Bsor = value) %>%
  tibble::rownames_to_column("ID")  # need the ID column to allow proper join (since it is a unique value)
bat.betadiversity.sim <- bat.simA1 %>%
  dplyr::select(row, col, regions, Bsim = value) %>%
  tibble::rownames_to_column("ID")
bat.betadiversity.sne <- bat.nesA1 %>%
  dplyr::select(row, col, regions, Bnes = value) %>%
  tibble::rownames_to_column("ID")

bat.betadiversity <- bat.betadiversity.sor %>%
  inner_join(bat.betadiversity.sim, by = "ID") %>%
  inner_join(bat.betadiversity.sne, by = "ID") %>%
  dplyr::select(regions, row, col, Bsor, Bsim, Bnes) #%>%
#filter(row != col) # this removes the self-self comparisons (e.g. AGO_Okavango - AGO_Okavango) which will obviously have betadiversity = 0

#------------------------------------------------------------------------------------------------------------------------------------------------------
#rats Bsor (total beta diversity)
rat.sorA <- melt(as.matrix(rat.beta$beta.sor), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to rat.sorA (above)
ratsite.reduced1 <- ratsite %>%
  dplyr::select(Location, row.Biogeogr = Biogeogr) # using select to also rename biogeogr column
ratsite.reduced2 <- ratsite %>%
  dplyr::select(Location, col.Biogeogr = Biogeogr) # using select to also rename biogeogr column

# now join the biogeographic regions to the betadiversity data.frame rat.sorA
rat.sorA1 <- rat.sorA %>%
  left_join(ratsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(ratsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

rat.sor.output <- rat.sorA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
rat.sor.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#rats Bsim (turnover beta diversity)
rat.simA <- melt(as.matrix(rat.beta$beta.sim), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to rat.simA (done above in rat.sor)
# now join the biogeographic regions to the betadiversity data.frame rat.simA
rat.simA1 <- rat.simA %>%
  left_join(ratsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(ratsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

rat.sim.output <- rat.simA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
rat.sim.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#rats Bnes (nestedness beta diversity)
rat.nesA <- melt(as.matrix(rat.beta$beta.sne), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to rat.sorA (above)
# now join the biogeographic regions to the betadiversity data.frame rat.nesA
rat.nesA1 <- rat.nesA %>%
  left_join(ratsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(ratsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

rat.nes.output <- rat.nesA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
rat.nes.output %>% print(n = Inf)  

# now create a data.frame with all three betadiversity values in a single data.frame and save to csv for further analysis and graphing
rat.betadiversity.sor <- rat.sorA1 %>%
  dplyr::select(row, col, regions, Bsor = value) %>%
  tibble::rownames_to_column("ID")  # need the ID column to allow proper join (since it is a unique value)
rat.betadiversity.sim <- rat.simA1 %>%
  dplyr::select(row, col, regions, Bsim = value) %>%
  tibble::rownames_to_column("ID")
rat.betadiversity.sne <- rat.nesA1 %>%
  dplyr::select(row, col, regions, Bnes = value) %>%
  tibble::rownames_to_column("ID")

rat.betadiversity <- rat.betadiversity.sor %>%
  inner_join(rat.betadiversity.sim, by = "ID") %>%
  inner_join(rat.betadiversity.sne, by = "ID") %>%
  dplyr::select(regions, row, col, Bsor, Bsim, Bnes) #%>%
  #filter(row != col) # this removes the self-self comparisons (e.g. AGO_Okavango - AGO_Okavango) which will obviously have betadiversity = 0

#------------------------------------------------------------------------------------------------------------------------------------------------------
#shrews Bsor (total beta diversity)
shrew.sorA <- melt(as.matrix(shrew.beta$beta.sor), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to shrew.sorA (above)
shrewsite.reduced1 <- shrewsite %>%
  dplyr::select(Location, row.Biogeogr = Biogeogr) # using select to also rename biogeogr column
shrewsite.reduced2 <- shrewsite %>%
  dplyr::select(Location, col.Biogeogr = Biogeogr) # using select to also rename biogeogr column

# now join the biogeographic regions to the betadiversity data.frame shrew.sorA
shrew.sorA1 <- shrew.sorA %>%
  left_join(shrewsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(shrewsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

shrew.sor.output <- shrew.sorA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
shrew.sor.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#shrews Bsim (turnover beta diversity)
shrew.simA <- melt(as.matrix(shrew.beta$beta.sim), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to shrew.simA (done above in shrew.sor)
# now join the biogeographic regions to the betadiversity data.frame shrew.simA
shrew.simA1 <- shrew.simA %>%
  left_join(shrewsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(shrewsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

shrew.sim.output <- shrew.simA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
shrew.sim.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#shrews Bnes (nestedness beta diversity)
shrew.nesA <- melt(as.matrix(shrew.beta$beta.sne), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to shrew.sorA (above)
# now join the biogeographic regions to the betadiversity data.frame shrew.nesA
shrew.nesA1 <- shrew.nesA %>%
  left_join(shrewsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(shrewsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

shrew.nes.output <- shrew.nesA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
shrew.nes.output %>% print(n = Inf)  

# now create a data.frame with all three betadiversity values in a single data.frame and save to csv for further analysis and graphing
shrew.betadiversity.sor <- shrew.sorA1 %>%
  dplyr::select(row, col, regions, Bsor = value) %>%
  tibble::rownames_to_column("ID")  # need the ID column to allow proper join (since it is a unique value)
shrew.betadiversity.sim <- shrew.simA1 %>%
  dplyr::select(row, col, regions, Bsim = value) %>%
  tibble::rownames_to_column("ID")
shrew.betadiversity.sne <- shrew.nesA1 %>%
  dplyr::select(row, col, regions, Bnes = value) %>%
  tibble::rownames_to_column("ID")

shrew.betadiversity <- shrew.betadiversity.sor %>%
  inner_join(shrew.betadiversity.sim, by = "ID") %>%
  inner_join(shrew.betadiversity.sne, by = "ID") %>%
  dplyr::select(regions, row, col, Bsor, Bsim, Bnes) #%>%
#filter(row != col) # this removes the self-self comparisons (e.g. AGO_Okavango - AGO_Okavango) which will obviously have betadiversity = 0

#------------------------------------------------------------------------------------------------------------------------------------------------------
# save all taxonomic beta diversity output as csv files

# this combines all components (Bsor, Bsim, Bnes) into a single csv file, separately for each group
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/R analysis/Community analyses/R Output/Betadiversity all_FINAL")

write.csv(bat.betadiversity, "bat.taxonomic.betadiversity.csv")
write.csv(rat.betadiversity, "rat.taxonomic.betadiversity.csv")
write.csv(shrew.betadiversity, "shrew.taxonomic.betadiversity.csv")

# this saves each component separately for each group
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/R analysis/Community analyses/R Output/Taxonomic betadiversity - between sites across all regions")

write.csv(bat.sor.output, "bat.sor.between all sites.csv")
write.csv(bat.sim.output, "bat.sim.between all sites.csv")
write.csv(bat.nes.output, "bat.nes.between all sites.csv")

write.csv(rat.sor.output, "rat.sor.between all sites.csv")
write.csv(rat.sim.output, "rat.sim.between all sites.csv")
write.csv(rat.nes.output, "rat.nes.between all sites.csv")

write.csv(shrew.sor.output, "shrew.sor.between all sites.csv")
write.csv(shrew.sim.output, "shrew.sim.between all sites.csv")
write.csv(shrew.nes.output, "shrew.nes.between all sites.csv")

#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------
# NOW REPEAT THIS FOR FUNCTIONAL BETA DIVERSITY

# First prepare traits data by reducing to four dimensions as required by betapart functional.beta.pair()

# bats
bat.pca <- prcomp(bats.traits[3:11], scale = TRUE)  # This excludes habitat features. Use [3:15] to include all traits
summary(bat.pca)

# convert the sites-species matrix into dataframe and extract only unique values for species - to get species list
bats.dataframe <- as.data.frame(bat_pres_abs) %>%    # converts matrix to dataframe
  tibble::rownames_to_column("Location") %>%   # add a new variable called Location based on rownames
  gather(key = Species, value = pres.abs, - Location) %>% # does a "reverse pivot" changing data to "long" format
  distinct(Species) # unique values only

# now extract pca coordinates for the first three axes  
bat.pca$x
bat.pca.df <- as.data.frame(bat.pca$x)  # convert to dataframe
bat.pca.df$Species <- bats.traits[ ,2]  # add species from bats.traits
bat.trait.tab1 <- bats.dataframe %>%  # this is now ready to be entered into betapart functional.beta.pair() as "traits"
  left_join(bat.pca.df, by = "Species") %>%
  dplyr::select(PC1, PC2, PC3, Species)

# now make species rownames
bat.trait.tab <- bat.trait.tab1 %>% column_to_rownames("Species")
bat.trait.table <- as.matrix(bat.trait.tab)

# Now prepare the corresponding bat records. But need to remove sites with less than 4 species (betapart requires that species richness be greater than the number of traits used; in this case 3)

# this is a useful way to filter out (ie remove) more than one item
'%ni%' <- Negate('%in%') # this creates the function %ni%, which works the opposite of %in% (which allows filtering of multiple items)

# need to remove ZAR_Kalahari (note the ZAR error) and ETH_Simien, both of which had 3 species
bat_pres_abs1 <- as.data.frame(bat_pres_abs) %>% tibble::rownames_to_column("Region") %>% 
  filter(Region %ni% c("ZAR_Kalahari", "ETH_Simien")) %>% 
  tibble::column_to_rownames("Region")
bat_pres_abs2 <- as.matrix(bat_pres_abs1)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Now repeat this for rats
rat.pca <- prcomp(rats.traits[3:10], scale = TRUE)  # This excludes habitat features. Use [3:15] to include all traits
summary(bat.pca)

# convert the sites-species matrix into dataframe and extract only unique values for species - to get species list
rats.dataframe <- as.data.frame(rat_pres_abs) %>%    # converts matrix to dataframe
  tibble::rownames_to_column("Location") %>%   # add a new variable called Location based on rownames
  gather(key = Species, value = pres.abs, - Location) %>% # does a "reverse pivot" changing data to "long" format
  distinct(Species) # unique values only

# now extract pca coordinates for the first three axes  
rat.pca$x
rat.pca.df <- as.data.frame(rat.pca$x)  # convert to dataframe
rat.pca.df$Species <- rats.traits[ ,2]  # add species from rats.traits
rat.trait.tab1 <- rats.dataframe %>%  # this is now ready to be entered into betapart functional.beta.pair() as "traits"
  left_join(rat.pca.df, by = "Species") %>%
  dplyr::select(PC1, PC2, PC3, Species)

# now make species rownames
rat.trait.tab <- rat.trait.tab1 %>% column_to_rownames("Species")
rat.trait.table <- as.matrix(rat.trait.tab)

# Now prepare the corresponding rat records. But need to remove sites with less than 4 species (betapart requires that species richness be greater than the number of traits used; in this case 3)

# this is a useful way to filter out (ie remove) more than one item
'%ni%' <- Negate('%in%') # this creates the function %ni%, which works the opposite of %in% (which allows filtering of multiple items)

# need to remove 18 sites had 3 species or less
rat_pres_abs1 <- as.data.frame(rat_pres_abs) %>% tibble::rownames_to_column("Region") %>% 
  filter(Region %ni% c("CMR_Kupe2", "CMR_Kupe3", "ETH_Bale3", "ETH_Bale4", "ETH_Simiens1", "KEN_Turkana2", "LBR_Nimba1", "MLI_Mopti", 
                       "NAM_Gobabeb", "NAM_Gorrasis", "NAM_Tumasberg", "TZA_KilimanjaroA4", "TZA_KilimanjaroA5", "ZAF_Karoo1", "ZAF_Karoo2", 
                       "ZAF_Karoo3", "ZAF_Karoo4", "ZAF_Karoo5")) %>% 
  tibble::column_to_rownames("Region")
rat_pres_abs2 <- as.matrix(rat_pres_abs1)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# And now shrews
shrew.pca <- prcomp(shrews.traits[3:10], scale = TRUE)  # This excludes habitat features. Use [3:15] to include all traits
summary(bat.pca)

# convert the sites-species matrix into dataframe and extract only unique values for species - to get species list
shrews.dataframe <- as.data.frame(shrew_pres_abs) %>%    # converts matrix to dataframe
  tibble::rownames_to_column("Location") %>%   # add a new variable called Location based on rownames
  gather(key = Species, value = pres.abs, - Location) %>% # does a "reverse pivot" changing data to "long" format
  distinct(Species) # unique values only

# now extract pca coordinates for the first three axes  
shrew.pca$x
shrew.pca.df <- as.data.frame(shrew.pca$x)  # convert to dataframe
shrew.pca.df$Species <- shrews.traits[ ,2]  # add species from rats.traits
shrew.trait.tab1 <- shrews.dataframe %>%  # this is now ready to be entered into betapart functional.beta.pair() as "traits"
  left_join(shrew.pca.df, by = "Species") %>%
  dplyr::select(PC1, PC2, Species)  # note only using 2 axes because of the low species richness per site

# now make species rownames
shrew.trait.tab <- shrew.trait.tab1 %>% column_to_rownames("Species")
shrew.trait.table <- as.matrix(shrew.trait.tab)

# Now prepare the corresponding shrew records. But need to remove sites with less than 2 species (betapart requires that species richness be greater than the number of traits used; in this case 1)

# this is a useful way to filter out (ie remove) more than one item
'%ni%' <- Negate('%in%') # this creates the function %ni%, which works the opposite of %in% (which allows filtering of multiple items)

# for 2 PCA axes - need to remove 67 sites had 2 species or less
shrew_pres_abs1 <- as.data.frame(shrew_pres_abs) %>% tibble::rownames_to_column("Region") %>%
  filter(Region %ni% c("CMR_Kupe1", "ETH_Alatish", "ETH_Bale1", "ETH_Bale3", "ETH_Bale4", "ETH_Bale5", "ETH_Chebera", "ETH_Hugumburda",
                       "ETH_Hunkolo", "ETH_Jiren", "ETH_Kaka", "ETH_Simiens2", "ETH_Simiens3", "ETH_Simiens4", "ETH_Yerer", "KEN_Kasigau",
                       "KEN_Meru", "KEN_Mukogodo", "Ken_Nairobi1", "KEN_Shimba", "KEN_Turkana1", "KEN_Turkana2", "LBR_Nimba1", "MLI_Mopti",
                       "MOZ_Cabora", "MOZ_Gorongosa", "MRT_Chott", "MWI_Lengwe", "MWI_Liwonde", "MWI_Zomba", "NAM_Keetmanshoop", "NAM_Tumasberg",
                       "NGA_Gambari", "SEN_Sabodala", "SLE_Bumbuna1", "SLE_Bumbuna2", "SOM_Jubba", "SWZ_Edwaleni", "SWZ_Maguga", "SWZ_Mlawula",
                       "TCD_Zakouma", "TZA_Kilimanjaro3", "TZA_Kingu", "UGA_Rwenzori4", "ZAF_Asante", "ZAF_Baviaanskloof", "ZAF_Karkloof",
                       "ZAF_Karoo2", "ZAF_Karoo3", "ZAF_Karoo4", "ZAF_Karoo5", "ZAF_Karoo6", "ZAF_Knersvlakte", "ZAF_Mankwe", "ZAF_Mariepskop",
                       "ZAF_MtZebra", "ZAF_Nossob", "ZAF_Nyslvley", "ZAF_Pretoriuskop", "ZAF_Rolfontein", "ZAF_Satara", "ZAF_Sneeuberg",
                       "ZAF_Swartberg", "ZAF_Wakefield", "ZAF_WillemPretorius", "ZMB_Kafue1", "ZWE_Sengwa")) %>%
  tibble::column_to_rownames("Region")

# for 1 PCA axis - need to remove ?? sites had 1 species (less than 2 species) - but this doesn't work; perhaps needs at least 2 traits?
# shrew_pres_abs1 <- as.data.frame(shrew_pres_abs) %>% tibble::rownames_to_column("Region") %>% 
#   filter(Region %ni% c("CMR_Kupe1", "ETH_Bale1", "ETH_Bale3", "ETH_Bale4", "ETH_Bale5", "ETH_Chebera", "ETH_Hugumburda", 
#                        "ETH_Jiren", "ETH_Simiens2", "ETH_Simiens3", "ETH_Simiens4", "ETH_Yerer", "KEN_Meru", "Ken_Nairobi1", "KEN_Turkana1", 
#                        "KEN_Turkana2", "LBR_Nimba1", "MLI_Mopti", "MOZ_Cabora", "MOZ_Gorongosa", "MWI_Lengwe", "MWI_Liwonde", "MWI_Zomba", 
#                        "NAM_Keetmanshoop", "NAM_Tumasberg", "SLE_Bumbuna1", "SLE_Bumbuna2", "SOM_Jubba", "SWZ_Edwaleni", "SWZ_Maguga", "SWZ_Mlawula",
#                        "TCD_Zakouma", "TZA_Kilimanjaro3", "TZA_Kingu", "ZAF_Asante", "ZAF_Baviaanskloof", "ZAF_Karkloof", "ZAF_Karoo2", 
#                        "ZAF_Karoo3", "ZAF_Karoo4", "ZAF_Karoo5", "ZAF_MtZebra", "ZAF_Nossob", "ZAF_Nyslvley", "ZAF_Satara", "ZAF_Swartberg", 
#                        "ZAF_WillemPretorius", "ZMB_Kafue1")) %>% 
#   tibble::column_to_rownames("Region")

shrew_pres_abs2 <- as.matrix(shrew_pres_abs1)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Functional beta diversity calculations using betapart beta.pair(), which compares between all sites (i.e. does not combine sites into a single output value)

bat.beta.func <- functional.beta.pair(bat_pres_abs2, bat.trait.table, index.family = "sorensen")
rat.beta.func <- functional.beta.pair(rat_pres_abs2, rat.trait.table, index.family = "sorensen")
shrew.beta.func <- functional.beta.pair(shrew_pres_abs2, shrew.trait.table, index.family = "sorensen")

#------------------------------------------------------------------------------------------------------------------------------------------------------
#bats Bsor (total functional beta diversity)
bat.fsorA <- melt(as.matrix(bats.func.Bsor$funct.beta.sor), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# now join the biogeographic regions to the functional betadiversity data.frame bat.fsorA
bat.fsorA1 <- bat.fsorA %>%
  left_join(batsite.reduced1, by = c('row' = 'Location')) %>%  # batsite.reduced1 and 2 created way above
  left_join(batsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

bat.fsor.output <- bat.fsorA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
bat.fsor.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#bats Bsim (turnover functional beta diversity)
bat.fsimA <- melt(as.matrix(bats.func.Bsor$funct.beta.sim), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to bat.fsimA (done above in bat.sor)
# now join the biogeographic regions to the betadiversity data.frame bat.fsimA
bat.fsimA1 <- bat.fsimA %>%
  left_join(batsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(batsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

bat.fsim.output <- bat.fsimA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
bat.fsim.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#bats Bnes (nestedness functional beta diversity)
bat.fnesA <- melt(as.matrix(bats.func.Bsor$funct.beta.sne), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to bat.sorA (above)
# now join the biogeographic regions to the betadiversity data.frame bat.nesA
bat.fnesA1 <- bat.fnesA %>%
  left_join(batsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(batsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

bat.fnes.output <- bat.fnesA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
bat.fnes.output %>% print(n = Inf)  

# now create a data.frame with all three betadiversity values in a single data.frame and save to csv for further analysis and graphing
bat.fbetadiversity.sor <- bat.fsorA1 %>%
  dplyr::select(row, col, regions, Bsor = value) %>%
  tibble::rownames_to_column("ID")  # need the ID column to allow proper join (since it is a unique value)
bat.fbetadiversity.sim <- bat.fsimA1 %>%
  dplyr::select(row, col, regions, Bsim = value) %>%
  tibble::rownames_to_column("ID")
bat.fbetadiversity.sne <- bat.fnesA1 %>%
  dplyr::select(row, col, regions, Bnes = value) %>%
  tibble::rownames_to_column("ID")

bat.fbetadiversity <- bat.fbetadiversity.sor %>%
  inner_join(bat.fbetadiversity.sim, by = "ID") %>%
  inner_join(bat.fbetadiversity.sne, by = "ID") %>%
  dplyr::select(regions, row, col, Bsor, Bsim, Bnes) #%>%
#filter(row != col) # this removes the self-self comparisons (e.g. AGO_Okavango - AGO_Okavango) which will obviously have betadiversity = 0

#------------------------------------------------------------------------------------------------------------------------------------------------------
#rats Bsor (total functional beta diversity)
rat.fsorA <- melt(as.matrix(rats.func.Bsor$funct.beta.sor), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# now join the biogeographic regions to the functional betadiversity data.frame rat.fsorA
rat.fsorA1 <- rat.fsorA %>%
  left_join(ratsite.reduced1, by = c('row' = 'Location')) %>%  # ratsite.reduced1 and 2 created way above
  left_join(ratsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

rat.fsor.output <- rat.fsorA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
rat.fsor.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#rats Bsim (turnover functional beta diversity)
rat.fsimA <- melt(as.matrix(rats.func.Bsor$funct.beta.sim), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to rat.fsimA (done above in rat.sor)
# now join the biogeographic regions to the betadiversity data.frame rat.fsimA
rat.fsimA1 <- rat.fsimA %>%
  left_join(ratsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(ratsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

rat.fsim.output <- rat.fsimA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
rat.fsim.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#rats Bnes (nestedness functional beta diversity)
rat.fnesA <- melt(as.matrix(rats.func.Bsor$funct.beta.sne), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to rat.sorA (above)
# now join the biogeographic regions to the betadiversity data.frame rat.nesA
rat.fnesA1 <- rat.fnesA %>%
  left_join(ratsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(ratsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

rat.fnes.output <- rat.fnesA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
rat.fnes.output %>% print(n = Inf)  

# now create a data.frame with all three betadiversity values in a single data.frame and save to csv for further analysis and graphing
rat.fbetadiversity.sor <- rat.fsorA1 %>%
  dplyr::select(row, col, regions, Bsor = value) %>%
  tibble::rownames_to_column("ID")  # need the ID column to allow proper join (since it is a unique value)
rat.fbetadiversity.sim <- rat.fsimA1 %>%
  dplyr::select(row, col, regions, Bsim = value) %>%
  tibble::rownames_to_column("ID")
rat.fbetadiversity.sne <- rat.fnesA1 %>%
  dplyr::select(row, col, regions, Bnes = value) %>%
  tibble::rownames_to_column("ID")

rat.fbetadiversity <- rat.fbetadiversity.sor %>%
  inner_join(rat.fbetadiversity.sim, by = "ID") %>%
  inner_join(rat.fbetadiversity.sne, by = "ID") %>%
  dplyr::select(regions, row, col, Bsor, Bsim, Bnes) #%>%
#filter(row != col) # this removes the self-self comparisons (e.g. AGO_Okavango - AGO_Okavango) which will obviously have betadiversity = 0

#------------------------------------------------------------------------------------------------------------------------------------------------------
#shrews Bsor (total functional beta diversity)
shrew.fsorA <- melt(as.matrix(shrews.func.Bsor$funct.beta.sor), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# now join the biogeographic regions to the functional betadiversity data.frame shrew.fsorA
shrew.fsorA1 <- shrew.fsorA %>%
  left_join(shrewsite.reduced1, by = c('row' = 'Location')) %>%  # shrewsite.reduced1 and 2 created way above
  left_join(shrewsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

shrew.fsor.output <- shrew.fsorA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
shrew.fsor.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#shrews Bsim (turnover functional beta diversity)
shrew.fsimA <- melt(as.matrix(shrews.func.Bsor$funct.beta.sim), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to shrew.fsimA (done above in shrew.sor)
# now join the biogeographic regions to the betadiversity data.frame shrew.fsimA
shrew.fsimA1 <- shrew.fsimA %>%
  left_join(shrewsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(shrewsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

shrew.fsim.output <- shrew.fsimA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
shrew.fsim.output %>% print(n = Inf)  # this allows you to see all the rows of an object. or use a number instead of Inf to define number of rows

#shrews Bnes (nestedness functional beta diversity)
shrew.fnesA <- melt(as.matrix(shrews.func.Bsor$funct.beta.sne), varnames = c("row", "col")) # takes the betadiversity matrix and converts into data.frame

# create the biogeographic regions for row and column separately to then join to shrew.sorA (above)
# now join the biogeographic regions to the betadiversity data.frame shrew.nesA
shrew.fnesA1 <- shrew.fnesA %>%
  left_join(shrewsite.reduced1, by = c('row' = 'Location')) %>%
  left_join(shrewsite.reduced2, by = c('col' = 'Location')) %>%
  unite("regions", row.Biogeogr:col.Biogeogr, sep="-") #%>%   # this combines the two columns into one using unite()

shrew.fnes.output <- shrew.fnesA1 %>%
  group_by(regions) %>%
  summarise(mean = mean(value), SD = sd(value))
shrew.fnes.output %>% print(n = Inf)  

# now create a data.frame with all three betadiversity values in a single data.frame and save to csv for further analysis and graphing
shrew.fbetadiversity.sor <- shrew.fsorA1 %>%
  dplyr::select(row, col, regions, Bsor = value) %>%
  tibble::rownames_to_column("ID")  # need the ID column to allow proper join (since it is a unique value)
shrew.fbetadiversity.sim <- shrew.fsimA1 %>%
  dplyr::select(row, col, regions, Bsim = value) %>%
  tibble::rownames_to_column("ID")
shrew.fbetadiversity.sne <- shrew.fnesA1 %>%
  dplyr::select(row, col, regions, Bnes = value) %>%
  tibble::rownames_to_column("ID")

shrew.fbetadiversity <- shrew.fbetadiversity.sor %>%
  inner_join(shrew.fbetadiversity.sim, by = "ID") %>%
  inner_join(shrew.fbetadiversity.sne, by = "ID") %>%
  dplyr::select(regions, row, col, Bsor, Bsim, Bnes) #%>%
#filter(row != col) # this removes the self-self comparisons (e.g. AGO_Okavango - AGO_Okavango) which will obviously have betadiversity = 0

#------------------------------------------------------------------------------------------------------------------------------------------------------
# save all functional beta diversity output as csv files

# this combines all components (Bsor, Bsim, Bnes) into a single csv file, separately for each group
setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/R analysis/Community analyses/R Output/Betadiversity all_FINAL")

write.csv(bat.fbetadiversity, "bat.functional.betadiversity.csv")
write.csv(rat.fbetadiversity, "rat.functional.betadiversity.csv")
write.csv(shrew.fbetadiversity, "shrew.functional.betadiversity.csv")

# this saves each component separately for each group

setwd("C:/Dropbox/Ara/Publications/Mss - current/Bats - Turnover comparisons with rodents/R analysis/Community analyses/R Output/Functional betadiversity - between sites across all regions")

# for bats
write.csv(bat.fsor.output, "bat.func_sor.between all sites.csv")
write.csv(bat.fsim.output, "bat.func_sim.between all sites.csv")
write.csv(bat.fnes.output, "bat.func_nes.between all sites.csv")

# for rats
write.csv(rat.fsor.output, "rat.func_sor.between all sites.csv")
write.csv(rat.fsim.output, "rat.func_sim.between all sites.csv")
write.csv(rat.fnes.output, "rat.func_nes.between all sites.csv")

# for shrews
write.csv(shrew.fsor.output, "shrew.func_sor.between all sites.csv")
write.csv(shrew.fsim.output, "shrew.func_sim.between all sites.csv")
write.csv(shrew.fnes.output, "shrew.func_nes.between all sites.csv")
