###################################################################
library("SingleCellExperiment")
library("iSEE")
library("shiny")
###########################################
sce_small <- readRDS("small_sce.rds")
#class(sce_small)
#sce_small@colData
iSEE(sce_small)
###########################################

