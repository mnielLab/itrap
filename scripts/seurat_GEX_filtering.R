#!/usr/bin/env Rscript

#################
### LIBRAREIS ###
#################

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(hdf5r)
#library(data.table)
#library(R.utils)

################################
###  COMMAND LINE ARGUMENTS  ###
################################

# RUN
#/home/tuba-nobackup/shared/R/R-3.6.1/bin/Rscript script.R input.rds output.csv

cArgs <- commandArgs(TRUE)

input_file   <- cArgs[1] # 'raw_feature_bc_matrix.h5'
out_file     <- cArgs[2] # 'gex_filtering.txt'
violin_plot  <- cArgs[3] #
scatter_plot <- cArgs[4] #


# Check if the correct number of input arguments were given
if (length(cArgs) < 4 | length(cArgs) > 4) {
  stop("The number of input arguments is not correct. There should be 2 input arguments! 
       Input arguments are: input.rds and output.csv")
}

#################
###    LOAD   ###
#################
gex.data <- Read10X_h5(input_file, use.names = TRUE, unique.features = TRUE)

#################
###    MAIN   ###
#################
# Initialize the Seurat object with the raw (non-normalized data).
#gex <- CreateSeuratObject(counts = gex.data[['Gene Expression']], project = "gex",
#                          min.cells = 0,
#                          min.features = 10) # min.features=200 'Antibody Capture'

gex <- tryCatch({CreateSeuratObject(counts = gex.data[['Gene Expression']], project = "gex",
                                    min.cells = 0,min.features = 10)},
                error = function(e){CreateSeuratObject(counts = gex.data, project = "gex",
                                                       min.cells = 0,min.features = 10)})

gex[["percent.mt"]] <- PercentageFeatureSet(gex, pattern = "^MT-")

#################
###    PLOT   ###
#################
# Visualize QC metrics as a violin plot
VlnPlot(gex, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(file=violin_plot, device='png', dpi=100)


apply(gex@assays$RNA@counts,2,function(x)(100*max(x))/sum(x)) -> gex$Percent.Largest.Gene

# Write QC metrics for follow-up plotting
as_tibble(gex[[c("nCount_RNA","nFeature_RNA","percent.mt","Percent.Largest.Gene")]],rownames="Cell.Barcode") -> qc.metrics

min_feat = 200
max_feat = 2500
max_mito = 10
min_rna = 1000
max_rna = 12000

# Plot UMI count versus gene count per GEM. Color by mitochondrial load
qc.metrics %>%
  arrange(percent.mt) %>%
  ggplot(aes(nCount_RNA,nFeature_RNA,colour=percent.mt)) + 
  geom_point() + 
  scale_color_gradientn(colors=c("black","blue","green2","yellow","red")) +
  ggtitle("QC metrics") +
  geom_hline(yintercept = min_feat) +
  geom_hline(yintercept = max_feat)
ggsave(file=scatter_plot, device='png', dpi=100)
      
###################
###    FILTER   ###
###################

#gex <- subset(gex, subset = nCount_RNA > min_rna & nCount_RNA < max_rna & nFeature_RNA > min_features)
#nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5
#gex <- subset(gex, subset = nCount_RNA > 1000 & nCount_RNA < 12000 & percent.mt < 10 & Percent.Largest.Gene < 15)
gex <- subset(gex, subset = nFeature_RNA > min_feat & nFeature_RNA < max_feat & percent.mt < max_mito)


write(colnames(gex),file=out_file)
      
      
# Write data to output
#data_to_write_out <- as.data.frame(as.matrix(data@meta.data)) #@scale.data
#write.csv(data_to_write_out, out_file, quote=FALSE)
