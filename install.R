if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("iSEE", dependencies = TRUE)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("SingleCellExperiment")

install.packages("shiny")

install.packages("tidyverse")
install.packages("rmarkdown")
install.packages("httr")
install.packages("shinydashboard")
install.packages('leaflet')
