#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library("SingleCellExperiment")
library("iSEE")
library("shiny")
library(ggplot2)
#sce_small <- readRDS("sce_updatedSeuObj_all_classes_curated.RDS");
#sce_small <- load("SCE_DLPFC-n3_tran-etal.rda");
# Define UI for application that draws a histogram
library(scRNAseq)
sce <- ReprocessedAllenData(assays = "tophat_counts")   # specifying the assays to speed up the example
#sce
ui <- fluidPage(
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

                
                #sce_small <- readRDS("/Users/mandala2/OneDrive - National Institutes of Health/shinyR_v1/sce_updatedSeuObj_all_classes_curated.RDS")
                
                #class(sce_small);
                #sce_small@colData
                iSEE(sce)
                
                #library(rsconnect)
                #rsconnect::deployApp('/Users/mandala2/OneDrive - National Institutes of Health/shinyR_v1/')
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$distPlot <- renderPlot({
        # generate bins based on input$bins from ui.R
        x    <- faithful[, 2]
        bins <- seq(min(x), max(x), length.out = input$bins + 1)

        # draw the histogram with the specified number of bins
        hist(x, breaks = bins, col = 'darkgray', border = 'white')
    })
}

# Run the application 
shinyApp(ui = ui, server = server)

