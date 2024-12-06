library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(scLinear)
library(scMRMA)
library(tibble)
library(viridis)
set.seed(42)
options(future.globals.maxSize = 2 * 1024^3) # 2 GiB


# Get data ----------------------------------------------------------------

# Create Surat object using Load10X_Spatial
data_dir <- "./local/spatial/lung"
h5_file <- "CytAssist_11mm_FFPE_Human_Lung_Cancer_filtered_feature_bc_matrix.h5"
lung <- Load10X_Spatial(data_dir, filename = h5_file)

# Data preprocessing ------------------------------------------------------

# Check hetrogeniety
plot1 <- VlnPlot(lung, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(lung, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

# Normalize with SCTransform to preserve biological variance
lung <- SCTransform(lung, assay = "Spatial", verbose = TRUE)

saveRDS(lung,"./local/spatial/lung/lung_normalized.rds")

# Dimensionality reduction, clustering and visualization ---------------------

lung <- readRDS("./local/spatial/lung/lung_normalized.rds")

# Run PCA, find neighbors, clustering and UMAP
lung <- RunPCA(lung, assay = "SCT", verbose = FALSE)
lung <- FindNeighbors(lung, reduction = "pca", dims = 1:30)
lung <- FindClusters(lung, verbose = FALSE)
lung <- RunUMAP(lung, reduction = "pca", dims = 1:30)

# Visualize clustering results
p1 <- DimPlot(lung, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(lung, label = TRUE, label.size = 3)
p1 + p2
SpatialDimPlot(lung, cells.highlight = CellsByIdentities(
  object = lung, idents = c(2, 1, 4, 3, 5, 8)), 
  facet.highlight = TRUE, ncol = 3)

# Cell annotation using scMRMA
DefaultAssay(lung) <- "SCT"
load(system.file("data", "Human_PanglaoDB.Rdata", package = "scMRMA"))
annotation <- scMRMA(input = lung, species = "Hs")

# Add cell type to the Seurat object (annotation level can be changed)
lung@meta.data$cell_type <- annotation$multiR$annotationResult[["Level1"]]

# Plot UMAP and spatial plot containing the cell types
p1<- DimPlot(lung, group.by = "cell_type", reduction = "umap")
p2 <- SpatialDimPlot(lung, group.by = "cell_type", images = "slice1")
p1 + p2

saveRDS(lung,"./local/spatial/lung/lung_annotated.rds")

# Predict ADT values ------------------------------------------------------

lung <- readRDS("./local/spatial/lung/lung_annotated.rds")

# Load the trained model
pipe <- create_adt_predictor()
pipe <- load_pretrained_model(pipe, file = "./local/lung/trained_model.joblib")

# Predict ADT values
lung@assays[["predicted_ADT"]] <-  adt_predict(pipe = pipe,
                                             gexp = lung@assays[["Spatial"]],
                                             normalize = TRUE)

saveRDS(lung,"./local/spatial/lung/lung_predicted.rds")


# Plot predicted ADT values for the different cell types ----------------------

lung <- readRDS("./local/spatial/lung/lung_predicted.rds")

# Create data frame of predicted ADT values
adt_predicted_sclinear <- as.matrix(lung@assays[["predicted_ADT"]]@data)
adt_predicted_sclinear <- adt_predicted_sclinear %>% as.data.frame() %>% rownames_to_column("gene") %>% pivot_longer(!gene,names_to = "cell", values_to = "predicted_sclinear")

# Make a map to change the protein names
adt_map <- c(
  "CD3-CD3E" = "CD3", 
  "CD4-CD4" = "CD4", 
  "CD8a-CD8A" = "CD8", 
  "CD11c-ITGAX" = "CD11c", 
  "CD14-CD14" = "CD14", 
  "CD15-FUT4" = "CD15", 
  "CD16-FCGR3A" = "CD16", 
  "CD19-CD19" = "CD19", 
  "CD127-IL7R" = "CD127", 
  "CD25-IL2RA" = "CD25", 
  "CD56-NCAM1" = "CD56", 
  "CD45RO" = "CD45RO", 
  "CD45RA" = "CD45RA", 
  "Mouse-IgG2a" = "Mouse-IgG2a", 
  "Mouse-IgG1" = "Mouse-IgG1", 
  "Mouse-IgG2b" = "Mouse-IgG2b", 
  "PD-1-PDCD1" = "PD-1", 
  "TIGIT-TIGIT" = "TIGIT", 
  "CD171-L1CAM" = "CD171", 
  "HER2-ERBB2" = "HER-2", 
  "CD140a-PDGFRA" = "CD140a", 
  "CD138-SDC1" = "CD138", 
  "EGFR-EGFR" = "EGFR", 
  "CD326-EPCAM" = "CD326", 
  "PD-L1-CD274" = "PD-L1", 
  "Podoplanin-PDPN" = "PDPN", 
  "Podocalyxin-PODXL" = "PODXL", 
  "CD133-PROM1" = "CD133", 
  "CD324-CDH1" = "CD324", 
  "CD244-CD244" = "CD244", 
  "CD44-CD44" = "CD44", 
  "CD10-MME" = "CD10"
)

# Replace names of proteins
adt_predicted_sclinear$gene <- adt_map[adt_predicted_sclinear$gene]

# Get cell types for each cell
meta <- lung@meta.data %>% rownames_to_column("cell") %>% dplyr::select("cell", "cell_type")

# Create data frame of predicted ADT values
DF <- adt_predicted_sclinear %>% full_join(meta, by = c("cell"))
DF <- DF %>% arrange(gene)
DF$gene <- factor(DF$gene, levels = unique(DF$gene))

write.table(DF, file = "./local/spatial/lung/bit_table_lung.csv", sep = ",", col.names = TRUE, row.names = FALSE)

lung <- read.table("./local/spatial/lung/bit_table_lung.csv", header = T, sep=',')
lung <- lung %>% distinct()

# Calculate predicted ADT value for each protein per cell type and save as csv
for (type in unique(lung$cell_type)) {
  adt_predicted_summary <- lung %>% filter(cell_type == type) %>% group_by(gene) %>% dplyr::summarise(mean = mean(predicted_sclinear))
  write.table(adt_predicted_summary, file = paste0("./local/spatial/lung/", type, ".csv"), sep = ",", col.names = TRUE, row.names = FALSE)
}

# Create plots for each cell type
for (type in unique(lung$cell_type)) {
  # Read csv-file containing mean predicted ADT values for each protein
  adt_predicted_celltype <- read.table(paste0("./local/spatial/lung/", type, ".csv"), header = TRUE, sep = ",")
  
  # Plot the predicted ADT values for each protein per cell type
  ggplot(adt_predicted_celltype, aes(x=reorder(gene, mean), y= mean,fill=mean)) + 
    geom_bar(stat="identity", col="black") + 
    coord_flip() +
    theme_classic() + 
    scale_fill_gradientn(colours = inferno(11)) + 
    ylab("Mean predicted ADT value (normalized)") + 
    xlab("Protein") + 
    theme(legend.position = "none") +
    ggtitle(paste0("Lung Cancer Visium CytAssist \n", type)) + 
    theme(plot.title = element_text(hjust = 0.5)) +
    expand_limits(y = 1)
  
  ggsave(paste0("./local/spatial/lung/figures/lung_spatial_", type, ".pdf"), width = 20/2, height = 29/2, units = 'cm')
  ggsave(paste0("./local/spatial/lung/figures/lung_spatial_", type, ".png"), width = 20/2, height = 29/2, units = 'cm')
}


# Plot predicted ADT values for each protein - spatial plot ---------------

lung <- readRDS("./local/spatial/lung/lung_predicted.rds")

# Map to change the protein names
adt_map <- c(
  "CD3-CD3E" = "CD3", 
  "CD4-CD4" = "CD4", 
  "CD8a-CD8A" = "CD8", 
  "CD11c-ITGAX" = "CD11c", 
  "CD14-CD14" = "CD14", 
  "CD15-FUT4" = "CD15", 
  "CD16-FCGR3A" = "CD16", 
  "CD19-CD19" = "CD19", 
  "CD127-IL7R" = "CD127", 
  "CD25-IL2RA" = "CD25", 
  "CD56-NCAM1" = "CD56", 
  "CD45RO" = "CD45RO", 
  "CD45RA" = "CD45RA", 
  "Mouse-IgG2a" = "Mouse-IgG2a", 
  "Mouse-IgG1" = "Mouse-IgG1", 
  "Mouse-IgG2b" = "Mouse-IgG2b", 
  "PD-1-PDCD1" = "PD-1", 
  "TIGIT-TIGIT" = "TIGIT", 
  "CD171-L1CAM" = "CD171", 
  "HER2-ERBB2" = "HER-2", 
  "CD140a-PDGFRA" = "CD140a", 
  "CD138-SDC1" = "CD138", 
  "EGFR-EGFR" = "EGFR", 
  "CD326-EPCAM" = "CD326", 
  "PD-L1-CD274" = "PD-L1", 
  "Podoplanin-PDPN" = "PDPN", 
  "Podocalyxin-PODXL" = "PODXL", 
  "CD133-PROM1" = "CD133", 
  "CD324-CDH1" = "CD324", 
  "CD244-CD244" = "CD244", 
  "CD44-CD44" = "CD44", 
  "CD10-MME" = "CD10"
)

# Change from gene names to protein names
DefaultAssay(lung) <- "Spatial" # Set default assay to something else than predicted_ADT
adt_data <- lung@assays[["predicted_ADT"]]@data # Extract the data matrix from the predicted_ADT assay
gene_names <- rownames(adt_data) # Get current gene names
new_gene_names <- sapply(gene_names, function(g) adt_map[[g]] %||% g) # Map the gene names using adt_map
rownames(adt_data) <- new_gene_names # Update with protein names
lung[["predicted_ADT"]] <- NULL # Remove the old predicted_ADT assay
lung[["predicted_ADT"]] <- CreateAssayObject(data = adt_data, key = "adt_") # Add the new assay with correct names
saveRDS(lung, "./local/spatial/lung/lung_predicted_protein_names.rds")


# Read RDS with protein names
lung <- readRDS("./local/spatial/lung/lung_predicted_protein_names.rds")

# Set default assay to be predicted_ADT
DefaultAssay(lung) <- "predicted_ADT"

# Make list of the protein names
proteins <- unique(rownames(lung@assays[["predicted_ADT"]]@data))

# Plot predicted ADT values (spatial) for each protein
for (protein in proteins) {
  SpatialFeaturePlot(lung, features = c(protein))
  ggsave(paste0("./local/spatial/lung/figures/", protein, "_lung_spatial.pdf"), width = 20/2, height = 29/2, units = 'cm')
  ggsave(paste0("./local/spatial/lung/figures/", protein, "_lung_spatial.png"), width = 20/2, height = 29/2, units = 'cm')
}

# Plot predicted ADt values (spatial) for each protein together with cell annotation
for (protein in proteins) {
  p1 <- SpatialFeaturePlot(lung, features = c(protein))
  p2 <- SpatialDimPlot(lung, group.by = "cell_type", images = "slice1")
  p1 + p2
  ggsave(paste0("./local/spatial/lung/figures/annotation_", protein, "_lung_spatial.pdf"), width = 15, height = 10, units = 'cm')
  ggsave(paste0("./local/spatial/lung/figures/annotation_", protein, "_lung_spatial.png"), width = 15, height = 10, units = 'cm')
}


# Create dotplot ----------------------------------------------------------

# Read RDS with protein names
lung <- readRDS("./local/spatial/lung/lung_predicted_protein_names.rds")

# Set the default assay to the predicted_ADT
DefaultAssay(lung) <- "predicted_ADT"

# Create a DotPlot for selected proteins across cell types
selected_proteins <- c("CD326", "CD324", "EGFR", "CD14", "CD3", "CD19", "CD56")

DotPlot(lung, features = selected_proteins, group.by = "cell_type") +
  scale_color_gradientn(colours = inferno(11)) +
  theme_minimal() + 
  labs(title = "Predicted ADT Values", x = "Proteins", y = "Cell Type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the dot plot
ggsave("./local/spatial/lung/figures/dotplot_predicted_ADT.pdf", width = 15, height = 10, units = "cm")
ggsave("./local/spatial/lung/figures/dotplot_predicted_ADT.png", width = 15, height = 10, units = "cm")






