# Final beta diversity comparisons, using the output from betapart (saved as csv files), Ara Monadjem (14 May 2021)

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
library(pairwise.adonis) # not working properly (for me)
library(reshape2)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# input data community data. These are beta diversities calculated using betapart beta.pair (ie site by site, and not than region by region)

# taxonomic betadiversity (pair-wise comparisons by sites across all regions)

bat.tax <- read.csv("bat.taxonomic.betadiversity.csv")
rat.tax <- read.csv("rat.taxonomic.betadiversity.csv")
shrew.tax <- read.csv("shrew.taxonomic.betadiversity.csv")

bat.func <- read.csv("bat.functional.betadiversity.csv")
rat.func <- read.csv("rat.functional.betadiversity.csv")
shrew.func <- read.csv("shrew.functional.betadiversity.csv")

bat.phylo <- read.csv("bat.phylo.betadiversity.csv")
rat.phylo <- read.csv("rat.phylo.betadiversity.csv")
shrew.phylo <- read.csv("shrew.phylo.betadiversity.csv")

#------------------------------------------------------------------------------------------------------------------------------------------------------
# rework the data for taxonomic betadiversity

# this removes the self-self comparisons (e.g. AGO_Okavango - AGO_Okavango) which will obviously have betadiversity = 0
bat.tax1 <- bat.tax %>% filter(row != col) %>% 
  separate(regions, c("rowBio","colBio"), sep = "-") %>%  # need to separate regions to select unique combinations
  filter(!duplicated(paste0(pmax(rowBio, colBio), pmin(rowBio, colBio)))) %>% # this selects unique combinations ie. Congolian-Ethiopian and Ethiopian-Congolian
  unite(regions, rowBio:colBio, sep = "-") %>% # return the regions as before
  group_by(regions) %>% 
  summarise(Bsor = mean(Bsor), Bsim = mean(Bsim), Bnes = mean(Bnes)) %>% 
  mutate_if(is.numeric, round, digits=2)
#bat.tax1 %>% print(n = Inf)

rat.tax1 <- rat.tax %>% filter(row != col) %>% 
  separate(regions, c("rowBio","colBio"), sep = "-") %>%  # need to separate regions to select unique combinations
  filter(!duplicated(paste0(pmax(rowBio, colBio), pmin(rowBio, colBio)))) %>% # this selects unique combinations ie. Congolian-Ethiopian and Ethiopian-Congolian
  unite(regions, rowBio:colBio, sep = "-") %>% # return the regions as before
  group_by(regions) %>% 
  summarise(Bsor = mean(Bsor), Bsim = mean(Bsim), Bnes = mean(Bnes)) %>% 
  mutate_if(is.numeric, round, digits=2)
#rat.tax1 %>% print(n = Inf)

shrew.tax1 <- shrew.tax %>% filter(row != col) %>% 
  separate(regions, c("rowBio","colBio"), sep = "-") %>%  # need to separate regions to select unique combinations
  filter(!duplicated(paste0(pmax(rowBio, colBio), pmin(rowBio, colBio)))) %>% # this selects unique combinations ie. Congolian-Ethiopian and Ethiopian-Congolian
  unite(regions, rowBio:colBio, sep = "-") %>% # return the regions as before
  group_by(regions) %>% 
  summarise(Bsor = mean(Bsor), Bsim = mean(Bsim), Bnes = mean(Bnes)) %>% 
  mutate_if(is.numeric, round, digits=2)
#shrew.tax1 %>% print(n = Inf)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# rework the data for functional betadiversity

# this removes the self-self comparisons (e.g. AGO_Okavango - AGO_Okavango) which will obviously have betadiversity = 0
bat.func1 <- bat.func %>% filter(row != col) %>% 
  separate(regions, c("rowBio","colBio"), sep = "-") %>%  # need to separate regions to select unique combinations
  filter(!duplicated(paste0(pmax(rowBio, colBio), pmin(rowBio, colBio)))) %>% # this selects unique combinations ie. Congolian-Ethiopian and Ethiopian-Congolian
  unite(regions, rowBio:colBio, sep = "-") %>% # return the regions as before
  group_by(regions) %>% 
  summarise(Bsor = mean(Bsor), Bsim = mean(Bsim), Bnes = mean(Bnes)) %>% 
  mutate_if(is.numeric, round, digits=2)
#bat.func1 %>% print(n = Inf)

rat.func1 <- rat.func %>% filter(row != col) %>% 
  separate(regions, c("rowBio","colBio"), sep = "-") %>%  # need to separate regions to select unique combinations
  filter(!duplicated(paste0(pmax(rowBio, colBio), pmin(rowBio, colBio)))) %>% # this selects unique combinations ie. Congolian-Ethiopian and Ethiopian-Congolian
  unite(regions, rowBio:colBio, sep = "-") %>% # return the regions as before
  group_by(regions) %>% 
  summarise(Bsor = mean(Bsor), Bsim = mean(Bsim), Bnes = mean(Bnes)) %>% 
  mutate_if(is.numeric, round, digits=2)
#rat.func1 %>% print(n = Inf)

shrew.func1 <- shrew.func %>% filter(row != col) %>% 
  separate(regions, c("rowBio","colBio"), sep = "-") %>%  # need to separate regions to select unique combinations
  filter(!duplicated(paste0(pmax(rowBio, colBio), pmin(rowBio, colBio)))) %>% # this selects unique combinations ie. Congolian-Ethiopian and Ethiopian-Congolian
  unite(regions, rowBio:colBio, sep = "-") %>% # return the regions as before
  group_by(regions) %>% 
  summarise(Bsor = mean(Bsor), Bsim = mean(Bsim), Bnes = mean(Bnes)) %>% 
  mutate_if(is.numeric, round, digits=2)
#shrew.func1 %>% print(n = Inf)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# rework the data for phylogenetic betadiversity

# this removes the self-self comparisons (e.g. AGO_Okavango - AGO_Okavango) which will obviously have betadiversity = 0
bat.phylo1 <- bat.phylo %>% filter(row != col) %>% 
  separate(regions, c("rowBio","colBio"), sep = "-") %>%  # need to separate regions to select unique combinations
  filter(!duplicated(paste0(pmax(rowBio, colBio), pmin(rowBio, colBio)))) %>% # this selects unique combinations ie. Congolian-Ethiopian and Ethiopian-Congolian
  unite(regions, rowBio:colBio, sep = "-") %>% # return the regions as before
  group_by(regions) %>% 
  summarise(Bsor = mean(Bsor), Bsim = mean(Bsim), Bnes = mean(Bnes)) %>% 
  mutate_if(is.numeric, round, digits=2)
#bat.phylo1 %>% print(n = Inf)

rat.phylo1 <- rat.phylo %>% filter(row != col) %>% 
  separate(regions, c("rowBio","colBio"), sep = "-") %>%  # need to separate regions to select unique combinations
  filter(!duplicated(paste0(pmax(rowBio, colBio), pmin(rowBio, colBio)))) %>% # this selects unique combinations ie. Congolian-Ethiopian and Ethiopian-Congolian
  unite(regions, rowBio:colBio, sep = "-") %>% # return the regions as before
  group_by(regions) %>% 
  summarise(Bsor = mean(Bsor), Bsim = mean(Bsim), Bnes = mean(Bnes)) %>% 
  mutate_if(is.numeric, round, digits=2)
#rat.phylo1 %>% print(n = Inf)

shrew.phylo1 <- shrew.phylo %>% filter(row != col) %>% 
  separate(regions, c("rowBio","colBio"), sep = "-") %>%  # need to separate regions to select unique combinations
  filter(!duplicated(paste0(pmax(rowBio, colBio), pmin(rowBio, colBio)))) %>% # this selects unique combinations ie. Congolian-Ethiopian and Ethiopian-Congolian
  unite(regions, rowBio:colBio, sep = "-") %>% # return the regions as before
  group_by(regions) %>% 
  summarise(Bsor = mean(Bsor), Bsim = mean(Bsim), Bnes = mean(Bnes)) %>% 
  mutate_if(is.numeric, round, digits=2)
#shrew.phylo1 %>% print(n = Inf)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Combining data into a single dataframe object; which requires tidying data first and removing duplicate values (e.g. BDI_Bururi-AGO_Kavango and AGO_Kavango-BDI_Bururi)

# taxonomic betadiversity
bat.tax2 <- bat.tax %>% filter(row != col) %>% 
  filter(!duplicated(paste0(pmax(as.character(row), as.character(col)), pmin(as.character(row), as.character(col))))) %>%
  unite(sites, row:col, sep = "-") %>% pivot_longer(Bsor:Bnes, names_to = "Component", values_to = "Value") %>% 
  mutate(group = "bats", betadiversity = "taxonomic")

rat.tax2 <- rat.tax %>% filter(row != col) %>% 
  filter(!duplicated(paste0(pmax(as.character(row), as.character(col)), pmin(as.character(row), as.character(col))))) %>%
  unite(sites, row:col, sep = "-") %>% pivot_longer(Bsor:Bnes, names_to = "Component", values_to = "Value") %>% 
  mutate(group = "rats", betadiversity = "taxonomic")

shrew.tax2 <- shrew.tax %>% filter(row != col) %>% 
  filter(!duplicated(paste0(pmax(as.character(row), as.character(col)), pmin(as.character(row), as.character(col))))) %>%
  unite(sites, row:col, sep = "-") %>% pivot_longer(Bsor:Bnes, names_to = "Component", values_to = "Value") %>% 
  mutate(group = "shrews", betadiversity = "taxonomic")

taxonomic <- bind_rows(bat.tax2, rat.tax2, shrew.tax2)
#------------------------------------------------------------------------------------------------------------------------------------------------------
# functional betadiversity: this groups 

bat.func2 <- bat.func %>% filter(row != col) %>% 
  filter(!duplicated(paste0(pmax(as.character(row), as.character(col)), pmin(as.character(row), as.character(col))))) %>%
  unite(sites, row:col, sep = "-") %>% pivot_longer(Bsor:Bnes, names_to = "Component", values_to = "Value") %>% 
  mutate(group = "bats", betadiversity = "functional")

rat.func2 <- rat.func %>% filter(row != col) %>% 
  filter(!duplicated(paste0(pmax(as.character(row), as.character(col)), pmin(as.character(row), as.character(col))))) %>%
  unite(sites, row:col, sep = "-") %>% pivot_longer(Bsor:Bnes, names_to = "Component", values_to = "Value") %>% 
  mutate(group = "rats", betadiversity = "functional")

shrew.func2 <- shrew.func %>% filter(row != col) %>% 
  filter(!duplicated(paste0(pmax(as.character(row), as.character(col)), pmin(as.character(row), as.character(col))))) %>%
  unite(sites, row:col, sep = "-") %>% pivot_longer(Bsor:Bnes, names_to = "Component", values_to = "Value") %>% 
  mutate(group = "shrews", betadiversity = "functional")

functional <- bind_rows(bat.func2, rat.func2, shrew.func2)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# functional betadiversity: this groups 

bat.phylo2 <- bat.phylo %>% filter(row != col) %>% 
  filter(!duplicated(paste0(pmax(as.character(row), as.character(col)), pmin(as.character(row), as.character(col))))) %>%
  unite(sites, row:col, sep = "-") %>% pivot_longer(Bsor:Bnes, names_to = "Component", values_to = "Value") %>% 
  mutate(group = "bats", betadiversity = "phylogenetic")

rat.phylo2 <- rat.phylo %>% filter(row != col) %>% 
  filter(!duplicated(paste0(pmax(as.character(row), as.character(col)), pmin(as.character(row), as.character(col))))) %>%
  unite(sites, row:col, sep = "-") %>% pivot_longer(Bsor:Bnes, names_to = "Component", values_to = "Value") %>% 
  mutate(group = "rats", betadiversity = "phylogenetic")

shrew.phylo2 <- shrew.phylo %>% filter(row != col) %>% 
  filter(!duplicated(paste0(pmax(as.character(row), as.character(col)), pmin(as.character(row), as.character(col))))) %>%
  unite(sites, row:col, sep = "-") %>% pivot_longer(Bsor:Bnes, names_to = "Component", values_to = "Value") %>% 
  mutate(group = "shrews", betadiversity = "phylogenetic")

phylogenetic <- bind_rows(bat.phylo2, rat.phylo2, shrew.phylo2)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Complete FINAL summary dataframe

betadiversity<- bind_rows(taxonomic, functional, phylogenetic)

write.csv(betadiversity, "betadiversity_all.csv")




##### END #####
#------------------------------------------------------------------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------

# taxonomic betadiversity (pair-wise comparisons by sites across all regions)
bat.sor <- read.csv("bat.sor.between all sites.csv")
bat.sim <- read.csv("bat.sim.between all sites.csv")
bat.nes <- read.csv("bat.nes.between all sites.csv")

rat.sor <- read.csv("rat.sor.between all sites.csv")
rat.sim <- read.csv("rat.sim.between all sites.csv")
rat.nes <- read.csv("rat.nes.between all sites.csv")

shrew.sor <- read.csv("shrew.sor.between all sites.csv")
shrew.sim <- read.csv("shrew.sim.between all sites.csv")
shrew.nes <- read.csv("shrew.nes.between all sites.csv")

# functional betadiversity (pair-wise comparisons by sites across all regions)
bat.fsor <- read.csv("bat.func_sor.between all sites.csv")
bat.fsim <- read.csv("bat.func_sim.between all sites.csv")
bat.fnes <- read.csv("bat.func_nes.between all sites.csv")

rat.fsor <- read.csv("rat.func_sor.between all sites.csv")
rat.fsim <- read.csv("rat.func_sim.between all sites.csv")
rat.fnes <- read.csv("rat.func_nes.between all sites.csv")

shrew.fsor <- read.csv("shrew.func_sor.between all sites.csv")
shrew.fsim <- read.csv("shrew.func_sim.between all sites.csv")
shrew.fnes <- read.csv("shrew.func_nes.between all sites.csv")

# functional betadiversity (pair-wise comparisons by sites across all regions)
bat.psor <- read.csv("bat.phylo_sor.between all sites.csv")
bat.psim <- read.csv("bat.phylo_sim.between all sites.csv")
bat.fnes <- read.csv("bat.phylo_nes.between all sites.csv")

rat.psor <- read.csv("rat.phylo_sor.between all sites.csv")
rat.psim <- read.csv("rat.phylo_sim.between all sites.csv")
rat.pnes <- read.csv("rat.phylo_nes.between all sites.csv")

shrew.psor <- read.csv("shrew.phylo_sor.between all sites.csv")
shrew.psim <- read.csv("shrew.phylo_sim.between all sites.csv")
shrew.pnes <- read.csv("shrew.phylo_nes.between all sites.csv")

#------------------------------------------------------------------------------------------------------------------------------------------------------
# Viewing the betadiversity values as a matrix

# taxonomic betadiversity
# bats
bat.sor.matrix <- bat.sor %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(bat.sor.matrix)

bat.sim.matrix <- bat.sim %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(bat.sim.matrix)

bat.nes.matrix <- bat.nes %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(bat.nes.matrix)

# rats
rat.sor.matrix <- rat.sor %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(rat.sor.matrix)

rat.sim.matrix <- rat.sim %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(rat.sim.matrix)

rat.nes.matrix <- rat.nes %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(rat.nes.matrix)

# shrews
shrew.sor.matrix <- shrew.sor %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(shrew.sor.matrix)

shrew.sim.matrix <- shrew.sim %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(shrew.sim.matrix)

shrew.nes.matrix <- shrew.nes %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(shrew.nes.matrix)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# functional betadiversity

# bats
bat.fsor.matrix <- bat.fsor %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(bat.fsor.matrix)

bat.fsim.matrix <- bat.fsim %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(bat.fsim.matrix)

bat.fnes.matrix <- bat.fnes %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(bat.fnes.matrix)

# rats
rat.fsor.matrix <- rat.fsor %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(rat.fsor.matrix)

rat.fsim.matrix <- rat.fsim %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(rat.fsim.matrix)

rat.fnes.matrix <- rat.fnes %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(rat.fnes.matrix)

# shrews
shrew.fsor.matrix <- shrew.fsor %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(shrew.fsor.matrix)

shrew.fsim.matrix <- shrew.fsim %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(shrew.fsim.matrix)

shrew.fnes.matrix <- shrew.fnes %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(shrew.fnes.matrix)

#------------------------------------------------------------------------------------------------------------------------------------------------------
# phylogenetic betadiversity

# bats
bat.psor.matrix <- bat.psor %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(bat.psor.matrix)

bat.psim.matrix <- bat.psim %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(bat.psim.matrix)

bat.pnes.matrix <- bat.pnes %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(bat.pnes.matrix)

# rats
rat.psor.matrix <- rat.psor %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(rat.psor.matrix)

rat.psim.matrix <- rat.psim %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(rat.psim.matrix)

rat.pnes.matrix <- rat.pnes %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(rat.pnes.matrix)

# shrews
shrew.psor.matrix <- shrew.psor %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(shrew.psor.matrix)

shrew.psim.matrix <- shrew.psim %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(shrew.psim.matrix)

shrew.pnes.matrix <- shrew.pnes %>% separate(regions, c("region1", "region2"), sep = "-", remove = T) %>% 
  dplyr::select(region1, region2, mean) %>% 
  pivot_wider(names_from = "region1", values_from = "mean")
as.matrix(shrew.pnes.matrix)
