# Summary of all Moran's Eigenvector Mapping analyses for the revision of the J Biogeogr manuscript, Ara Monadjem (26 June 2022)
# Nearly all of the original analyses have been done in other R scripts, and are referred to by name here

#-------------------------------------------------------------------------------------------------------------------------------------------
rm(list = ls()) # removes everything currently held in the R memory
graphics.off() # closes all open graphics windows
# to clear the screen click "control L"
# set.seed(123)

#-------------------------------------------------------------------------------------------------------------------------------------------
## Libraries:

library(tidyverse)
library(vegan)

#-------------------------------------------------------------------------------------------------------------------------------------------
# call data
# selection of significant MEMs (for bats, rats, shrews for each of taxonomic, functional, and phylogenetic)
# calculated in three R scripts: 'Bat_MEM_AM.R', 'Rat_MEM_AM.R', 'Shrew_MEM_AM.R'
load("bat.MEM.sel.beta.tax.RData")
load("bat.MEM.sel.beta.func.RData")
load("bat.MEM.sel.beta.phyl.RData")
load("rat.MEM.sel.beta.tax.RData")
load("rat.MEM.sel.beta.func.RData")
load("rat.MEM.sel.beta.phyl.RData")
load("shrew.MEM.sel.beta.tax.RData")
load("shrew.MEM.sel.beta.func.RData")
load("shrew.MEM.sel.beta.phyl.RData")

bat.MEM.sel.beta.tax
bat.MEM.sel.beta.func
bat.MEM.sel.beta.phyl

rat.MEM.sel.beta.tax
rat.MEM.sel.beta.func
rat.MEM.sel.beta.phyl

shrew.MEM.sel.beta.tax
shrew.MEM.sel.beta.func
shrew.MEM.sel.beta.phyl

#-------------------------------------------------------------------------------------------------------------------------------------------
# RDA of bat community on MEMs (not very useful?)

bat.rda.mem.up <- rda(batsp ~ ., data = bat.MEM.sel.beta.tax)
bat.rda.mem.lw <- rda(batsp ~ 1, data = bat.MEM.sel.beta.tax)

#bat.mods <- ordistep(bat.rda.mem.lw, scope = formula(bat.rda.mem.up),
#                     trace = FALSE)
# need a 2-step selection process and hence ordiR2step() - see below
bat.mods1 <- ordiR2step(bat.rda.mem.lw, bat.rda.mem.up, trace = FALSE)
bat.mods1$anova

# now test for significant environmental variables based on bat tax MEMs
bat.rda.mem.up1 <- rda(bat.MEM.sel.beta.tax ~ ., data = env.bat.final)
bat.rda.mem.lw1 <- rda(bat.MEM.sel.beta.tax ~ 1., data = env.bat.final)

# bat.modsa <- ordistep(bat.rda.mem.lw1, scope = formula(bat.rda.mem.up1),
#                      trace = FALSE)
# need a 2-step selection process and hence ordiR2step() - see below
bat.mods1a <- ordiR2step(bat.rda.mem.lw1, bat.rda.mem.up1, trace = FALSE)
bat.mods1a$anova

# now test for significant environmental variables based on bat tax MEMs
rat.rda.mem.up1 <- rda(rat.MEM.sel.beta.tax ~ ., data = env.rat.final)
rat.rda.mem.lw1 <- rda(rat.MEM.sel.beta.tax ~ 1., data = env.rat.final)

# rat.modsa <- ordistep(rat.rda.mem.lw1, scope = formula(rat.rda.mem.up1),
#                       trace = FALSE)
# need a 2-step selection process and hence ordiR2step() - see below
rat.mods1a <- ordiR2step(rat.rda.mem.lw1, rat.rda.mem.up1, trace = FALSE)
rat.mods1a$anova

# now test for significant environmental variables based on shrew tax MEMs
shrew.rda.mem.up1 <- rda(shrew.MEM.sel.beta.tax ~ ., data = env.shrew.final)
shrew.rda.mem.lw1 <- rda(shrew.MEM.sel.beta.tax ~ 1., data = env.shrew.final)

# rat.modsa <- ordistep(rat.rda.mem.lw1, scope = formula(rat.rda.mem.up1), 
#                       trace = FALSE)
# need a 2-step selection process and hence ordiR2step() - see below
shrew.mods1a <- ordiR2step(shrew.rda.mem.lw1, shrew.rda.mem.up1, trace = FALSE)
shrew.mods1a$anova

