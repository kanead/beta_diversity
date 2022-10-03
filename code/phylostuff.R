## The packages
library(ape)
library(caper)
library(phytools)

## A random tree (for the demo)
my_tree <- rcoal(15)

## A random dataset of two variables (for the PGLS)
my_data <- as.data.frame(matrix(rnorm(30), 15, 2))
rownames(my_data) <- my_tree$tip.label

## For matching the names in the data and the tree they are various options
# Can't find the package from Luis Verde (maybe I dreamt it) but they might be some useful info/scripts here: http://www.italian-journal-of-mammalogy.it/Good-practices-for-sharing-analysis-ready-data-in-mammalogy-and-biodiversity-research,101564,0,2.html

## For the exact matchings (i.e. are X names exactly in the tree)
# geiger::name.check # For matching one dataset to one tree
# dispRity::clean.data # For matching multiple datasets to multiple trees

###############
## Measuring Faith's phylogenetic diversity
###############

## The function for measuring Faith's PD per groups
#' @phy the phylogeny
#' @focal a vector of names of the focal group
#' @rest optional, a vector to scale from, can be nothing (not scaling), all the tip labels in phy (full scaling) or the tip labels in other communities (e.g. for calculating the phylogenetic beta diversity)
group.beta.phy <- function(phy, focal = NULL, rest = NULL) {
    
    if(is.null(focal)) {
        ## Faith's PD for the whole tree
        results <- sum(phy$edge.length)
    } else {
        ## Create a sub-tree
        to_drop     <- !(phy$tip.label %in% focal)
        sub_tree    <- drop.tip(phy, tip = phy$tip.label[to_drop])
        results <- sum(sub_tree$edge.length)
    }

    if(!is.null(rest)) {
        ## Create a sub-tree for rest
        to_drop    <- !(phy$tip.label %in% rest)
        sub_tree   <- drop.tip(phy, tip = phy$tip.label[to_drop])
        rest_brlen <- sum(sub_tree$edge.length)

        ## Dividing the group Faith PD by the total
        results <- 2*(sum(phy$edge.length))/(results + rest_brlen)
    }
    return(results)
}

## The total branch length
group.beta.phy(my_tree)

## The focal group tip labels
focal_group <- c("t1", "t2", "t3", "t4", "t5")

## The branch length for that group
group.beta.phy(my_tree, focal = focal_group)

## The other groups
rest_group  <- my_tree$tip.label[!c(my_tree$tip.label %in% focal_group)]

## The beta phylo diversity
abs(1 - group.beta.phy(my_tree, focal = focal_group, rest = rest_group))

community <- matrix(nrow = 2, ncol = Ntip(my_tree))
rownames(community) <- c("focal", "rest")
colnames(community) <- my_tree$tip.label
community["focal", ] <- as.integer(my_tree$tip.label %in% focal_group)
community["rest", ]  <- as.integer(my_tree$tip.label %in% rest_group)

## This:
BLij <- sum(my_tree$edge.length)
BLi  <- sum(drop.tip(my_tree, tip = my_tree$tip.label[!(my_tree$tip.label %in% focal_group)])$edge.length)
BLj  <- sum(drop.tip(my_tree, tip = my_tree$tip.label[!(my_tree$tip.label %in% rest_group)])$edge.length)
PhyloSor <- (2*BLij)/(BLi + BLj)
abs(1 - PhyloSor)

## Is equivalent to this
phylo.beta.pair(community, my_tree)$phylo.beta.sor




###############
## Measuring Pagel's Lamda for some traits
###############
#' @phy the phylogeny
#' @data the data to evaluate (one variable with species names)
#' @group optional, a vector of names of the focal group
#' @method either "K" for Blomberg's K or "lambda" for Pagel's lambda
#' @relative optional, whether to give the raw phylo signal or the relative one
group.phylosig <- function(phy, data, method, group = NULL, relative = FALSE) {
    
    if(is.null(group)) {
        ## Phylo signal for the whole tree
        phylo_sig <- phylosig(phy, data, method = method)[[1]]

    } else {
        ## Create a sub-tree
        to_drop <- !(phy$tip.label %in% group)
        sub_tree <- drop.tip(phy, tip = phy$tip.label[to_drop])
        sub_data <- data[c(!(names(data) %in% phy$tip.label[to_drop]))]
        phylo_sig <- phylosig(sub_tree, sub_data, method = method)[[1]]
    }

    if(relative) {
        ## Dividing the group's signal by the total one
        phylo_sig <- phylo_sig/phylosig(phy, data, method = method)[[1]]
    }
    return(phylo_sig)
}

## Selecting the variable (with named elements)
my_variable <- my_data[,1]
names(my_variable) <- rownames(my_data)

## Blomberg's K for the first variable
group.phylosig(my_tree, my_variable, method = "K")
## Pagel's lambda for the first variable
group.phylosig(my_tree, my_variable, method = "lambda")

## Blomberg's K for the focal group
group.phylosig(my_tree, my_variable, method = "K", group = focal_group)
## Pagel's lambda for the focal group
group.phylosig(my_tree, my_variable, method = "lambda", group = focal_group)


## Measuring Pagel's Lambda for a PGLS between two traits

## One of the column has to be species names
my_data <- cbind(my_data, my_tree$tip.label)
colnames(my_data) <-c("Variable1", "Variable2", "Names")

## Setting up the data
comp_data <- comparative.data(phy = my_tree, data = my_data, 
                              names.col = Names, vcv = TRUE, 
                              na.omit = FALSE, warn.dropped = TRUE)

## Running the PGLS
model_pgls <- pgls(Variable1 ~ Variable2, data = comp_data, lambda = "ML")

## Extracting the phylogenetic signal
summary(model_pgls)$param["lambda"]
## And the lambda detailes
summary(model_pgls)$param.CI["lambda"]





##function for joining up lists of databases

multi_join <- function(list_of_loaded_data, join_func, ...){
  
  require("dplyr")
  
  output <- Reduce(function(x, y) {join_func(x, y, ...)}, list_of_loaded_data)
  
  return(output)
}






