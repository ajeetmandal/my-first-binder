###################################################################
# library(Seurat)
# clustObject <- readRDS(file = "/Volumes/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/reprocess_11022018/all_01142019.RDS")
# 
# clustObject <- readRDS(file = "/Volumes/HBCC_analysis/analysis/kim2/single_nucleus_project_2018/Rstudio/re")
# ## check that Seurat package is at least v3.0
# utils::packageVersion('Seurat') < 3 
# ## check version of Seurat object 
# clustObject@version < 3
# 
# ## UpdateSeuratObject: Update old Seurat object to accommodate new features
# clustObject_updated <- UpdateSeuratObject(clustObject)
# class(clustObject)
# sce <- as.SingleCellExperiment(clustObject_updated)
# sce@colData
# saveRDS(sce, "/Users/mandala2/OneDrive - National Institutes of Health/shinyR_v1/sce_updatedSeuObj_all_classes_curated.RDS")
###########################################
# TRY
#iSEE(sce)
###########################################
library("SingleCellExperiment")
library("iSEE")
library("shiny")
library(ggplot2)

sce_small <- readRDS("/Users/mandala2/OneDrive - National Institutes of Health/shinyR_v1/sce_updatedSeuObj_all_classes_curated.RDS")
class(sce_small)
sce_small@colData
iSEE(sce_small)
