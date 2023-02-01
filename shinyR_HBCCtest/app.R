library("SingleCellExperiment")
library("iSEE")
library("shiny")
library("ggplot2")

#sce_small <- readRDS("small_sce.rds")
sce_small <- load("small_sce_v1.rds")
sce_small <- sce.dlpfc.tran
#sce_small <- sce_small[,1:1000]
# # Save a single object to a file
#saveRDS(sce_small, "/Users/mandala2/OneDrive - National Institutes of Health/shinyR_v4/small_sce_v1.rds")
#sce_small <- readRDS("small_sce_v1.rds")

#colData(sce_small)[,"sizeFactors(sce_small)"] <- sizeFactors(sce_small)
colormap <- ExperimentColorMap()
colormap <- synchronizeAssays(colormap, sce_small)
all_coordinates <- list()
custom_data_fun <- NULL
custom_stat_fun <- NULL

################################################################################
## Feature assay plot 1
################################################################################

plot.data <- data.frame(Y=assay(sce_small, 2, withDimnames=FALSE)[7711,], row.names = colnames(sce_small))
plot.data$X <- colData(sce_small)[,"cellType"];
plot.data$ColorBy <- colData(sce_small)[,"sum"];
plot.data <- subset(plot.data, !is.na(X) & !is.na(Y));

# Saving data for transmission
all_coordinates[['featAssayPlot1']] <- plot.data

# Setting up plot coordinates
plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- iSEE::jitterViolinPoints(plot.data$X, plot.data$Y, 
                                                width=0.4, varwidth=FALSE, adjust=1,
                                                method='quasirandom', nbins=NULL);

# Creating the plot
ggplot() +
  geom_violin(aes(x = X, y = Y, group = GroupBy), alpha = 0.2, data=plot.data, scale = 'width', width = 0.8) +
  geom_point(aes(y = Y, color = ColorBy, x = jitteredX), alpha = 1, plot.data, size=1) +
  labs(x = "Cluster", y = "PF4 (logcounts)", color = "log10_total_counts", title = "PF4 vs Cluster") +
  coord_cartesian(ylim = range(plot.data$Y, na.rm=TRUE), expand = TRUE) +
  scale_color_gradientn(colors=colDataColorMap(colormap, "log10_total_counts", discrete=FALSE)(21), na.value='grey50', limits=range(plot.data$ColorBy, na.rm=TRUE)) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(legend.position = 'bottom', legend.text=element_text(size=9),
        legend.title=element_text(size=11), legend.box = 'vertical',
        axis.text.x = element_text(angle=90, size=10, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=12), title=element_text(size=12))

################################################################################
## Column data plot 1
################################################################################

plot.data <- data.frame(Y = colData(sce_small)[,"sum"], row.names=colnames(sce_small));
plot.data$X <- colData(sce_small)[,"cellType"];
plot.data$ColorBy <- colData(sce_small)[,"detected"];
plot.data$ColorBy <- factor(plot.data$ColorBy);
plot.data <- subset(plot.data, !is.na(X) & !is.na(Y));

# Saving data for transmission
all_coordinates[['colDataPlot1']] <- plot.data

# Setting up plot coordinates
plot.data$GroupBy <- plot.data$X;
set.seed(100);
plot.data$jitteredX <- iSEE::jitterViolinPoints(plot.data$X, plot.data$Y, 
                                                width=0.4, varwidth=FALSE, adjust=1,
                                                method='quasirandom', nbins=NULL);

# Creating the plot
ggplot() +
  geom_violin(aes(x = X, y = Y, group = GroupBy), alpha = 0.2, data=plot.data, scale = 'width', width = 0.8) +
  geom_point(aes(y = Y, color = ColorBy, x = jitteredX), alpha = 1, plot.data, size=1) +
  labs(x = "Cluster", y = "log10_total_counts", color = "detected", title = "log10_total_counts vs Cluster") +
  coord_cartesian(ylim = range(plot.data$Y, na.rm=TRUE), expand = TRUE) +
  scale_color_manual(values=colDataColorMap(colormap, "detected", discrete=TRUE)(5198), na.value='grey50', drop=FALSE) +
  scale_fill_manual(values=colDataColorMap(colormap, "detected", discrete=TRUE)(3), na.value='grey50', drop=FALSE) +
  scale_x_discrete(drop = FALSE) +
  theme_bw() +
  theme(legend.position = 'bottom', legend.text=element_text(size=9),
        legend.title=element_text(size=11), legend.box = 'vertical',
        axis.text.x = element_text(angle=90, size=10, hjust=1, vjust=0.5),
        axis.text.y=element_text(size=10),
        axis.title=element_text(size=12), title=element_text(size=12))

################################################################################
## Heat map 1
################################################################################

value.mat <- as.matrix(assay(sce_small, 2)[1L, , drop=FALSE]);
plot.data <- reshape2::melt(value.mat, varnames = c('Y', 'X'));

plot.data[['OrderBy1']] <- factor(colData(sce_small)[['cellType']][match(plot.data$X, rownames(colData(sce_small)))]);
plot.data <- dplyr::arrange(plot.data, OrderBy1);
plot.data$X <- factor(plot.data$X, levels = unique(plot.data$X));

# Centering and scaling
plot.data$value <- plot.data$value - ave(plot.data$value, plot.data$Y);

# Creating the heat map
p0 <- ggplot(plot.data, aes(x = X, y = Y)) +
  geom_raster(aes(fill = value)) +
  labs(x='', y='') +
  scale_fill_gradientn(colors=c('purple','black','yellow'),
                       values=c(0,0.5,1),
                       limits=c(-5,5), na.value='grey50') +
  scale_y_discrete(expand=c(0, 0)) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.line=element_blank());
heatlegend <- cowplot::get_legend(p0 + theme(legend.position='bottom'));

# Adding annotations
legends <- list()

p1 <- ggplot(plot.data, aes(x = X, y = 1)) +
  geom_raster(aes(fill = OrderBy1)) +
  labs(x='', y='') +
  scale_y_continuous(breaks=1, labels='cellType') +
  scale_fill_manual(values=colDataColorMap(colormap, 'cellType', discrete=TRUE)(19), na.value='grey50', drop=FALSE, name='cellType') +
  theme(axis.text.x=element_blank(), axis.ticks=element_blank(), axis.title.x=element_blank(),
        rect=element_blank(), line=element_blank(), axis.title.y=element_blank(),
        plot.margin = unit(c(0,0,-0.5,0), 'lines'));
legends[[1]] <- cowplot::get_legend(p1 + theme(legend.position='bottom', plot.margin = unit(c(0,0,0,0), 'lines')));


################################################################################
## Reduced dimension plot 1
################################################################################

red.dim <- reducedDim(sce_small, 2);
plot.data <- data.frame(X = red.dim[, 1], Y = red.dim[, 2], row.names=colnames(sce_small));
plot.data$ColorBy <- colData(sce_small)[,"detected"];
plot.data$ColorBy <- factor(plot.data$ColorBy);
plot.data <- subset(plot.data, !is.na(X) & !is.na(Y));

# Saving data for transmission
all_coordinates[['redDimPlot1']] <- plot.data

# Creating the plot
ggplot() +
  geom_point(aes(x = X, y = Y, color = ColorBy), alpha = 1, plot.data, size=1) +
  labs(x = "Dimension 1", y = "Dimension 2", color = "detected", title = "(2) TSNE") +
  coord_cartesian(xlim = range(plot.data$X, na.rm = TRUE),
                  ylim = range(plot.data$Y, na.rm = TRUE), expand = TRUE) +
  scale_color_manual(values=colDataColorMap(colormap, "detected", discrete=TRUE)(5198), na.value='grey50', drop=FALSE) +
  scale_fill_manual(values=colDataColorMap(colormap, "detected", discrete=TRUE)(3), na.value='grey50', drop=FALSE) +
  theme_bw() +
  theme(legend.position = 'bottom', legend.box = 'vertical', legend.text=element_text(size=9), legend.title=element_text(size=11),
        axis.text=element_text(size=10), axis.title=element_text(size=12), title=element_text(size=12))




# Laying out the grid
cowplot::plot_grid(
  cowplot::plot_grid(
    p1 + theme(legend.position='none'),
    p0 + theme(legend.position='none'),
    ncol=1, align='v', rel_heights=c(0.1, 1)),
  heatlegend, ncol=1, rel_heights=c(0.9, 0.1))

################################################################################
## To guarantee the reproducibility of your code, you should also
## record the output of sessionInfo()
sessionInfo()
