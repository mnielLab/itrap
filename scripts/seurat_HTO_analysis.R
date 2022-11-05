#!/usr/bin/env Rscript

#################
### LIBRAREIS ###
#################

library(Seurat)
library(ggplot2)
library(data.table)
library(R.utils)

################################
###  COMMAND LINE ARGUMENTS  ###
################################

# RUN
#/home/tuba-nobackup/shared/R/R-3.6.1/bin/Rscript script.R input.rds output.csv

cArgs <- commandArgs(TRUE)

input_file   <- cArgs[1]
faux_file    <- cArgs[2]
out_file     <- cArgs[3]
ridge_plot   <- cArgs[4] # "hto_ridgeplot.png"
violin_plot  <- cArgs[5] # "hto_violinplot.png"
heatmap      <- cArgs[6] # "hto_heatmap.png"


# Check if the correct number of input arguments were given
if (length(cArgs) < 6 | length(cArgs) > 6) {
  stop("The number of input arguments is not correct. There should be 2 input arguments! 
       Input arguments are: input.rds and output.csv")
}

#################
###    LOAD   ###
#################
rds <- readRDS(input_file)

#################
###    MAIN   ###
#################
# Setup Seurat object
data <- CreateSeuratObject(counts=rds, assay='HTO')
# Normalize HTO data, here we use centered log-ratio (CLR) transformation
data <- NormalizeData(data, assay = "HTO", normalization.method = "CLR")

# If you have a very large dataset we suggest using kfunc = 'clara'. This is a k-medoid
# clustering function for large applications You can also play with additional parameters (see
# documentation for HTODemux()) to adjust the threshold for classification Here we are using
# the default settings
data <- HTODemux(data, assay = "HTO", positive.quantile = 0.99)

# Visualize
# Group cells based on the max HTO signal
Idents(data) <- "HTO_maxID"
RidgePlot(data, assay = "HTO", features = rownames(data[["HTO"]])[1:10], ncol = 10)
ggsave(file=ridge_plot, device='png', width=14, height=4)

# Only HTOs
Idents(data) <- "HTO_classification.global"
VlnPlot(data, features = "nCount_HTO", pt.size = 0.1, log = TRUE)
ggsave(file=violin_plot, device='png')

# HTO heatmap
# To increase the efficiency of plotting, you can subsample cells using the num.cells argument
HTOHeatmap(data, assay = "HTO") #, ncells = 5000
ggsave(file=heatmap, device='png', width=14, height=4) #, dpi=300


# Write data to output
data_to_write_out <- as.data.frame(as.matrix(data@meta.data)) #@scale.data
write.csv(data_to_write_out, out_file, quote=FALSE)
#fwrite(x = data_to_write_out, row.names = TRUE, file = out_file)