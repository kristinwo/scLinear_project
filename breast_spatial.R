library(dplyr)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(grid)
library(hdf5r)
library(Matrix)
library(Metrics)
library(pals)
library(patchwork)
library(reshape)
library(scLinear)
library(scMRMA)
library(Seurat)
library(SeuratData)
library(svglite)
library(tibble)
library(tidyverse)
library(viridis)
set.seed(42)
options(future.globals.maxSize = 600 * 1024^2)

# Get data ----------------------------------------------------------------

# Read h5 file
data_dir <- "./local/spatial/breast/"
h5_file <- "CytAssist_FFPE_Protein_Expression_Human_Breast_Cancer_filtered_feature_bc_matrix.h5"
breast_data <- Read10X_h5(paste0(data_dir, h5_file))

# Create Seurat object from the data in the h5 file
assay <- "RNA"
breast <- CreateSeuratObject(counts = breast_data[["Gene Expression"]], assay = assay, min.cells = 0, min.features = 0)
breast[["ADT"]] <- CreateAssayObject(breast_data[["Antibody Capture"]][, colnames(x = breast)])

# Add the image to the Seurat object
image <- Read10X_Image(
  image.dir = file.path(data_dir, 'spatial'),
  filter.matrix = TRUE,
  image.name = "tissue_lowres_image.png")

image <- image[Cells(x = breast)]
DefaultAssay(breast = image) <- assay
breast[["slice1"]] <- image


# Data preprocessing ------------------------------------------------------

# Check hetrogeniety
plot1 <- VlnPlot(breast, features = "nCount_RNA", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(breast, features = "nCount_RNA", images="slice1") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

# Normalize RNA with SCTransform to preserve biological variance
breast <- SCTransform(breast, assay = "RNA", verbose = FALSE)

# Normalize ADT with CLR
breast <- NormalizeData(breast, assay = "ADT", normalization.method = "CLR", margin = 2)

# Save normalized Seurat object
saveRDS(breast, "./local/spatial/breast/breast_normalized.rds")


# Dimensionality reduction, clustering and visualization ------------------

breast <- readRDS("./local/spatial/breast/breast_normalized.rds")

DefaultAssay(breast) <- "SCT"

# Run PCA, find neighbors, clustering and UMAP
breast <- RunPCA(breast, assay = "SCT", verbose = FALSE)
breast <- FindNeighbors(breast, reduction = "pca", dims = 1:30)
breast <- FindClusters(breast, verbose = FALSE)
breast <- RunUMAP(breast, reduction = "pca", dims = 1:30)

# Visualize clustering results
p1 <- DimPlot(breast, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(breast, label = TRUE, label.size = 3)
p1 + p2
SpatialDimPlot(breast, cells.highlight = CellsByIdentities(
  object = breast, idents = c(2, 1, 4, 3, 5, 8)), 
  facet.highlight = TRUE, ncol = 3)

# Cell annotation using scMRMA
DefaultAssay(breast) <- "SCT"
load(system.file("data", "Human_PanglaoDB.Rdata", package = "scMRMA"))
annotation <- scMRMA(input = breast, species = "Hs")

# Add cell type to the Seurat object (annotation level can be changed)
breast@meta.data$cell_type <- annotation$multiR$annotationResult[["Level2"]]

# Plot UMAP and spatial plot containing the cell types
p1<- DimPlot(breast, group.by = "cell_type", reduction = "umap")
p2 <- SpatialDimPlot(breast, group.by = "cell_type", images = "slice1")
p1 + p2

saveRDS(breast,"./local/spatial/breast/breast_annotated.rds")


# Predict ADT values ------------------------------------------------------

breast <- readRDS("./local/spatial/breast/breast_annotated.rds")

# Load the trained model
pipe <- create_adt_predictor()
pipe <- load_pretrained_model(pipe, file = "./local/spatial/breast/trained_model.joblib")

# Create a training and a test set
set.seed(42)
indx <- sample(1:length(colnames(breast)), size = length(colnames(breast)), replace = FALSE)
breast_train <- breast[,indx[1:2000]]
breast_test <- breast[,indx[2001:length(colnames(breast))]]
saveRDS(breast_train ,"./local/spatial/breast/breast_train.rds")
saveRDS(breast_test ,"./local/spatial/breast/breast_test.rds")

# Load training and test set
breast_train <- readRDS("./local/spatial/breast/breast_train.rds")
breast_test <- readRDS("./local/spatial/breast/breast_test.rds")

# Create predictor
pipe <- create_adt_predictor()

# Train predictor
pipe <- fit_predictor(pipe = pipe,
                      gexp_train = breast_train@assays[["RNA"]],
                      adt_train = breast_train@assays[["ADT"]],
                      normalize_gex = TRUE,
                      normalize_adt = TRUE)

# Save the trained model
save_trained_model(pipe = pipe, file = "./local/spatial/breast/trained_model.joblib")

# Load the trained model
pipe <- create_adt_predictor()
pipe <- load_pretrained_model(pipe, file = "./local/spatial/breast/trained_model.joblib")

# Evaluate predictor
eval_res <- evaluate_predictor(pipe = pipe,
                               gexp_test = breast_test@assays[["RNA"]],
                               adt_test = breast_test@assays[["ADT"]],
                               normalize_gex = TRUE,
                               normalize_adt = TRUE)

print(eval_res)

# Save evaluation results as a csv
write.table(eval_res, "./local/spatial/breast/evaluation_results_all_breast.csv")

# Predict ADT values
breast_test@assays[["predicted_ADT"]] <-  adt_predict(pipe = pipe,
                                                      gexp = breast_test@assays[["RNA"]],
                                                      normalize = TRUE)

saveRDS(breast,"./local/spatial/breast/breast_predicted.rds")


# Train a new model for each cell type ------------------------------------

# Load train and test data
breast_train <- readRDS("./local/spatial/breast/breast_train.rds")
breast_test <- readRDS("./local/spatial/breast/breast_test.rds")

# Find the common cell types for both datasets
cell_types_train <- unique(breast_train@meta.data[["cell_type"]])
cell_types_test <- unique(breast_test@meta.data[["cell_type"]])
cell_types <- intersect(cell_types_train, cell_types_test)

# Create a data frame to store the evaluation results for each cell type
results <- data.frame(RMSE = numeric(),
                      Pearson = numeric(),
                      Spearman = numeric(),
                      cell_type = character(),
                      stringsAsFactors = FALSE)


for (cell_type in cell_types) {
  
  # Subset the train and and test data for the specific cell type
  breast_train_subset <- subset(breast_train, subset = breast_train@meta.data[["cell_type"]] == cell_type)
  breast_test_subset <- subset(breast_test, subset = breast_test@meta.data[["cell_type"]] == cell_type)
  
  # Create predictor
  pipe <- create_adt_predictor()
  
  # Train predictor
  pipe <- fit_predictor(pipe = pipe,
                        gexp_train = breast_train_subset@assays[["RNA"]],
                        adt_train = breast_train_subset@assays[["ADT"]],
                        normalize_gex = TRUE,
                        normalize_adt = TRUE)
  
  # Save the trained model
  filename <- paste0("./local/spatial/breast/trained_model_", cell_type, ".joblib")
  save_trained_model(pipe = pipe, file = filename)
  
  # Load the trained model
  pipe <- create_adt_predictor()
  pipe <- load_pretrained_model(pipe, file = filename)
  
  # Evaluate predictor
  eval_res <- evaluate_predictor(pipe = pipe,
                                 gexp_test = breast_test_subset@assays[["RNA"]],
                                 adt_test = breast_test_subset@assays[["ADT"]],
                                 normalize_gex = TRUE,
                                 normalize_adt = TRUE)
  
  print(cell_type)
  print(eval_res)
  
  
  
  # Add the predicted adt assay
  breast_test_subset[["predicted_ADT"]] <-  adt_predict(pipe = pipe,
                                                        gexp = breast_test_subset@assays[["RNA"]],
                                                        normalize = TRUE)
  
  # Add the evaluation results to the results data frame
  results <- rbind(results, data.frame(RMSE = eval_res$RMSE,
                                       Pearson = eval_res$Pearson,
                                       Spearman = eval_res$Spearman,
                                       cell_type = cell_type))
}

# Write the evaluation results to a csv file
write.table(results, file = "./local/spatial/breast/evaluation_results_celltype_breast.csv", 
            sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Create correlation figure for each protein ------------------------------

breast <- readRDS("./local/spatial/breast/breast_predicted.rds")

# Create data frame of real and predicted ADT values
adt_real <- as.matrix(breast@assays[["ADT"]]@data)
adt_real <- adt_real %>% as.data.frame() %>% rownames_to_column("gene") %>% pivot_longer(!gene,names_to = "cell", values_to = "real")
adt_predicted_sclinear <- as.matrix(breast@assays[["predicted_ADT"]]@data)
adt_predicted_sclinear <- adt_predicted_sclinear %>% as.data.frame() %>% rownames_to_column("gene") %>% pivot_longer(!gene,names_to = "cell", values_to = "predicted_sclinear")

# Make a maps to change the protein names
adt_map <- c(
  "CD163.1" = "CD163",
  "CR2.1" = "CR2",
  "PCNA.1" = "PCNA",
  "VIM.1" = "VIM",
  "KRT5.1" = "KRT5",
  "CD68.1" = "CD68",
  "CEACAM8.1" = "CD66b",
  "PTPRC.1" = "CD45.1",
  "HLA-DRA" = "HLA-DR",
  "PAX5.1" = "PAX-5",
  "SDC1.1" = "CD138",
  "PTPRC.2" = "CD45.2",
  "CD8A.1" = "CD8",
  "BCL2.1" = "BCL2",
  "mouse-IgG2a" = "mouse-IgG2a",
  "mouse-IgG1k" = "mouse-IgG1k",
  "mouse-IgG2bk" = "mouse-IgG2bk",
  "rat-IgG2a" = "rat-IgG2a",
  "CD19.1" = "CD19",
  "PDCD1.1" = "PD-1",
  "ACTA2.1" = "ACTA2",
  "FCGR3A.1" = "CD16",
  "ITGAX.1" = "CD11c",
  "CXCR5.1" = "CD185",
  "EPCAM.1" = "CD326",
  "MS4A1.1" = "CD20",
  "CD3E.1" = "CD3",
  "CD14.1" = "CD14",
  "CD40.1" = "CD40",
  "PECAM1.1" = "CD31",
  "CD4.1" = "CD4",
  "ITGAM.1" = "CD11b",
  "CD27.1" = "CD27",
  "CCR7.1" = "CCR7",
  "CD274.1" = "PD-L1"
)


# Replace names of proteins
adt_real$gene <- adt_map[adt_real$gene]
adt_predicted_sclinear$gene <- adt_map[adt_predicted_sclinear$gene]

# Get cell types for each cell
meta <- breast@meta.data %>% rownames_to_column("cell") %>% dplyr::select("cell", "cell_type")

# Create data frame of real and predicted ADT values
DF <- adt_real %>% full_join(adt_predicted_sclinear, by = c("gene", "cell")) %>% full_join(meta, by = c("cell"))
DF <- DF %>% arrange(gene)
DF$gene <- factor(DF$gene, levels = unique(DF$gene))

write.table(DF, file = "./local/spatial/breast/bit_table_breast.csv", sep = ",", col.names = TRUE, row.names = FALSE)

breast <- read.table("./local/spatial/breast/bit_table_breast.csv", header = T, sep=',')
breast <- breast %>% distinct()

# Calculate pearson correlation for each protein
cors_pearson <- breast %>% group_by(gene) %>% dplyr::summarise(correlation = cor(real, predicted_sclinear, method="pearson"))
write.csv(cors_pearson, "./local/spatial/breast/pearson_correlation_proteins_breast.csv")

# Calculate spearman correlation for each protein
cors_spearman <- breast %>% group_by(gene) %>% dplyr::summarise(correlation = cor(real, predicted_sclinear, method="spearman"))
write.csv(cors_spearman, "./local/spatial/breast/spearman_correlation_proteins_breast.csv")

# Calculate RMSE for each protein
rmse <- breast %>% group_by(gene) %>% dplyr::summarise(rmse = rmse(real, predicted_sclinear))
write.csv(cors_spearman, "./local/spatial/breast/rmse_proteins_breast.csv")


# Plot the pearson correlation for each protein
ggplot(cors_pearson, aes(x=reorder(gene, correlation), y= correlation,fill=correlation)) + 
  geom_bar(stat="identity", col="black") + 
  coord_flip() +
  theme_classic() + 
  scale_fill_gradientn(colours = inferno(11)) + 
  ylab("Pearson\n(Real vs ScLinear)") + 
  xlab("Protein") + 
  theme(legend.position = "none") +
  ggtitle("Breast Cancer Visium CystAssist") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y = 1)

ggsave("./local/spatial/breast/figures/breast_Pearson_proteins.pdf", width = 20/2, height = 29/2, units = 'cm')
ggsave("./local/spatial/breast/figures/breast_Pearson_proteins.png", width = 20/2, height = 29/2, units = 'cm')

# Plot the spearman correlation for each protein
ggplot(cors_spearman, aes(x=reorder(gene, correlation), y= correlation,fill=correlation)) + 
  geom_bar(stat="identity", col="black") + 
  coord_flip() +
  theme_classic() + 
  scale_fill_gradientn(colours = inferno(11)) + 
  ylab("Spearman\n(Real vs ScLinear)") + 
  xlab("Protein") + 
  theme(legend.position = "none") +
  ggtitle("Breast Cancer Visium CystAssist") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y = 1)

ggsave("./local/spatial/breast/figures/breast_Spearman_proteins.pdf", width = 20/2, height = 29/2, units = 'cm')
ggsave("./local/spatial/breast/figures/breast_Spearman_proteins.png", width = 20/2, height = 29/2, units = 'cm')

# Plot the RMSE for each protein
ggplot(rmse, aes(x=reorder(gene, rmse), y= rmse,fill=rmse)) + 
  geom_bar(stat="identity", col="black") + 
  coord_flip() +
  theme_classic() + 
  scale_fill_gradientn(colours = inferno(11)) + 
  ylab("RMSE\n(Real vs ScLinear)") + 
  xlab("Protein") + 
  theme(legend.position = "none") +
  ggtitle("Breast Cancer Visium CystAssist") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y = 1)

ggsave("./local/spatial/breast/figures/breast_RMSE_proteins.pdf", width = 20/2, height = 29/2, units = 'cm')
ggsave("./local/spatial/breast/figures/breast_RMSE_proteins.png", width = 20/2, height = 29/2, units = 'cm')

# Summary of RMSE, Pearson and Spearman for breast spatial ----------------------

# Read csv files containing RMSE, Pearson and Spearman correlations
results_breast <- read.table("./local/spatial/breast/evaluation_results_all_breast.csv")

# Add a column specifying the tissue the data is from
results_breast$Tissue <- "Breast"

# Merge into one single data frame
results <- results_breast

# Reshape the data to long format for easier plotting
long_results <- melt(results, id.vars = "Tissue", 
                     measure.vars = c("RMSE", "Pearson", "Spearman"))

# Define colors for each tissue
colors <- c(
  "RMSE" = "#7D7D61", 
  "Pearson" = "#7B90AA",
  "Spearman" = "#BF878A"
)


# Create plot with RMSE, Pearson and Spearman correlation for spatial breast cancer
breast_plot <- ggplot(long_results, aes(x = variable, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.4) +
  labs(title = "Breast Cancer \nVisium CystAssist",
       y = "") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        legend.position = "none",
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 1)) + 
  scale_fill_manual(values = colors)
breast_plot

ggsave("./local/spatial/breast/figures/summary_statistics.pdf", width = 10, height = 8, units = 'cm')
ggsave("./local/spatial/breast/figures/summary_statistics.png", width = 10, height = 8, units = 'cm')


# RMSE, Pearson and Spearman for each cell type ---------------------------

# Read csv files containing RMSE, Pearson, and Spearman correlations
results_breast <- read.table("./local/spatial/breast/evaluation_results_all_breast.csv", header = TRUE)

# Read cell type data
results_breast_celltype <- read.csv("./local/spatial/breast/evaluation_results_celltype_breast.csv")

# Add a column to distinguish between overall and cell types
results_breast$cell_type <- "Overall"

# Combine overall and cell type data for each tissue
results_breast_combined <- rbind(results_breast, results_breast_celltype)

# Reshape data to long format for plotting
long_results_breast <- melt(results_breast_combined, id.vars = "cell_type", 
                            measure.vars = c("RMSE", "Pearson", "Spearman"))

# Define colors for each cell type in breast
colors_breast <- c(
  "B" = "#7B90AA",
  "Duct epithelial cell" = "#7D7D61",
  "Epithelial cell" = "#84A8A3",
  "Hematopoietic cell" = "#B1BECD",
  "Overall" = "#BF878A",
  "T" = "#C4B893"
)

# Plot cell types in breast cancer spatial
ggplot(long_results_breast, aes(x = variable, y = value, fill = cell_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.4) +
  labs(title = "Breast cancer \nVisium CyastAssist", y = "", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(0, 1)) + 
  scale_fill_manual(name = "Cell Type", values = colors_breast)

ggsave("./local/spatial/breast/figures/breast_celltypes.pdf", width = 15, height = 8, units = 'cm')
ggsave("./local/spatial/breast/figures/breast_celltypes.png", width = 15, height = 8, units = 'cm')

# Stacked barplot of cell count in each tissue ----------------------------

# Create level 1 annotation breast spatial data
breast <- readRDS("./local/spatial/breast/breast_normalized.rds")

DefaultAssay(breast) <- "SCT"

# Run PCA, find neighbors, clustering and UMAP
breast <- RunPCA(breast, assay = "SCT", verbose = FALSE)
breast <- FindNeighbors(breast, reduction = "pca", dims = 1:30)
breast <- FindClusters(breast, verbose = FALSE)
breast <- RunUMAP(breast, reduction = "pca", dims = 1:30)

# Visualize clustering results
p1 <- DimPlot(breast, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(breast, label = TRUE, label.size = 3)
p1 + p2
SpatialDimPlot(breast, cells.highlight = CellsByIdentities(
  object = breast, idents = c(2, 1, 4, 3, 5, 8)), 
  facet.highlight = TRUE, ncol = 3)

# Cell annotation using scMRMA
DefaultAssay(breast) <- "SCT"
load(system.file("data", "Human_PanglaoDB.Rdata", package = "scMRMA"))
annotation <- scMRMA(input = breast, species = "Hs")

# Add cell type to the Seurat object (annotation level can be changed)
breast@meta.data$cell_type <- annotation$multiR$annotationResult[["Level1"]]

# Plot UMAP and spatial plot containing the cell types
p1<- DimPlot(breast, group.by = "cell_type", reduction = "umap")
p2 <- SpatialDimPlot(breast, group.by = "cell_type", images = "slice1")
p1 + p2

saveRDS(breast,"./local/spatial/breast/breast_annotated_anno1.rds")

# Read level 1 annotated files
breast_sc <- readRDS("./local/breast/breast_prepared_anno1.rds")
breast_spatial <- readRDS("./local/spatial/breast/breast_annotated_anno1.rds")

# Extract cell type information and tissue labels
get_celltype_data <- function(seurat_obj, tissue) {
  data.frame(
    cell_type = seurat_obj@meta.data$cell_type,
    tissue = tissue
  )
}

# Combine data from all tissues
breast_sc_data <- get_celltype_data(breast_sc, "Breast CITE-seq")
breast_spatial_data <- get_celltype_data(breast_spatial, "Breast Visium CystAssist")
celltype_data <- bind_rows(breast_sc_data, breast_spatial_data)

# Calculate percentage of each cell type per tissue
celltype_percentages <- celltype_data %>%
  group_by(tissue, cell_type) %>% 
  dplyr::summarise(count = n()) %>% 
  group_by(tissue) %>%
  mutate(percentage = count / sum(count) * 100)

# Define colors for each cell type
colors <- c(
  "Hematopoietic cell" = "#7D7D61", 
  "Epithelial cell" = "#7B90AA",
  "Connective tissue cell" = "#BF878A"
)

# Plot with custom colors
ggplot(celltype_percentages, aes(x = tissue, y = percentage, fill = cell_type)) +
  geom_bar(stat = "identity") +
  labs(x = "", y = "Percentage of all cells") +
  scale_fill_manual(name = "Cell Type", values = colors) +
  theme_minimal() +
  theme(legend.position = "right")

ggsave("./local/spatial/breast/figures/cell_count_stacked_barplot.pdf", width = 11, height = 10, units = 'cm')
ggsave("./local/spatial/breast/figures/cell_count_stacked_barplot.png", width = 11, height = 10, units = 'cm')

# Plot ADT values for each protein - spatial plot ---------------

breast <- readRDS("./local/spatial/breast/breast_predicted.rds")

# Make a maps to change the protein names
adt_map <- c(
  "CD163.1" = "CD163",
  "CR2.1" = "CR2",
  "PCNA.1" = "PCNA",
  "VIM.1" = "VIM",
  "KRT5.1" = "KRT5",
  "CD68.1" = "CD68",
  "CEACAM8.1" = "CD66b",
  "PTPRC.1" = "CD45.1",
  "HLA-DRA" = "HLA-DR",
  "PAX5.1" = "PAX-5",
  "SDC1.1" = "CD138",
  "PTPRC.2" = "CD45.2",
  "CD8A.1" = "CD8",
  "BCL2.1" = "BCL2",
  "mouse-IgG2a" = "mouse-IgG2a",
  "mouse-IgG1k" = "mouse-IgG1k",
  "mouse-IgG2bk" = "mouse-IgG2bk",
  "rat-IgG2a" = "rat-IgG2a",
  "CD19.1" = "CD19",
  "PDCD1.1" = "PD-1",
  "ACTA2.1" = "ACTA2",
  "FCGR3A.1" = "CD16",
  "ITGAX.1" = "CD11c",
  "CXCR5.1" = "CD185",
  "EPCAM.1" = "CD326",
  "MS4A1.1" = "CD20",
  "CD3E.1" = "CD3",
  "CD14.1" = "CD14",
  "CD40.1" = "CD40",
  "PECAM1.1" = "CD31",
  "CD4.1" = "CD4",
  "ITGAM.1" = "CD11b",
  "CD27.1" = "CD27",
  "CCR7.1" = "CCR7",
  "CD274.1" = "PD-L1"
)

# Change from gene names to protein names
DefaultAssay(breast) <- "SCT" # Set default assay to something else than ADT
adt_data <- breast@assays[["ADT"]]@data # Extract the data matrix from the ADT assay
gene_names <- rownames(adt_data) # Get current gene names
new_gene_names <- sapply(gene_names, function(g) adt_map[[g]] %||% g) # Map the gene names using adt_map
rownames(adt_data) <- new_gene_names # Update with protein names
breast@assays[["ADT"]] <- NULL # Remove the old ADT assay
breast@assays[["ADT"]] <- CreateAssayObject(data = adt_data, key = "adt_") # Add the new assay with correct names
saveRDS(breast, "./local/spatial/breast/breast_predicted_protein_names_ADT.rds")


# Read RDS with protein names
breast <- readRDS("./local/spatial/breast/breast_predicted_protein_names_ADT.rds")

# Set default assay to be ADT
DefaultAssay(breast) <- "ADT"

# Make list of the protein names
proteins <- unique(rownames(breast@assays[["ADT"]]@data))

# Plot predicted ADT values (spatial) for each protein
for (protein in proteins) {
  SpatialFeaturePlot(breast, features = c(protein))
  ggsave(paste0("./local/spatial/breast/figures/ADT/", protein, "_breast_spatial.pdf"), width = 20/2, height = 29/2, units = 'cm')
  ggsave(paste0("./local/spatial/breast/figures/ADT/", protein, "_breast_spatial.png"), width = 20/2, height = 29/2, units = 'cm')
}

# Plot predicted ADt values (spatial) for each protein together with cell annotation
for (protein in proteins) {
  p1 <- SpatialFeaturePlot(breast, features = c(protein))
  p2 <- SpatialDimPlot(breast, group.by = "cell_type", images = "slice1")
  p1 + p2
  ggsave(paste0("./local/spatial/breast/figures/ADT/annotation_", protein, "_breast_spatial.pdf"), width = 15, height = 10, units = 'cm')
  ggsave(paste0("./local/spatial/breast/figures/ADT/annotation_", protein, "_breast_spatial.png"), width = 15, height = 10, units = 'cm')
}

# Plot predicted ADT values for each protein - spatial plot ---------------

breast <- readRDS("./local/spatial/breast/breast_predicted.rds")

# Make a maps to change the protein names
adt_map <- c(
  "CD163.1" = "CD163",
  "CR2.1" = "CR2",
  "PCNA.1" = "PCNA",
  "VIM.1" = "VIM",
  "KRT5.1" = "KRT5",
  "CD68.1" = "CD68",
  "CEACAM8.1" = "CD66b",
  "PTPRC.1" = "CD45.1",
  "HLA-DRA" = "HLA-DR",
  "PAX5.1" = "PAX-5",
  "SDC1.1" = "CD138",
  "PTPRC.2" = "CD45.2",
  "CD8A.1" = "CD8",
  "BCL2.1" = "BCL2",
  "mouse-IgG2a" = "mouse-IgG2a",
  "mouse-IgG1k" = "mouse-IgG1k",
  "mouse-IgG2bk" = "mouse-IgG2bk",
  "rat-IgG2a" = "rat-IgG2a",
  "CD19.1" = "CD19",
  "PDCD1.1" = "PD-1",
  "ACTA2.1" = "ACTA2",
  "FCGR3A.1" = "CD16",
  "ITGAX.1" = "CD11c",
  "CXCR5.1" = "CD185",
  "EPCAM.1" = "CD326",
  "MS4A1.1" = "CD20",
  "CD3E.1" = "CD3",
  "CD14.1" = "CD14",
  "CD40.1" = "CD40",
  "PECAM1.1" = "CD31",
  "CD4.1" = "CD4",
  "ITGAM.1" = "CD11b",
  "CD27.1" = "CD27",
  "CCR7.1" = "CCR7",
  "CD274.1" = "PD-L1"
)

# Change from gene names to protein names
DefaultAssay(breast) <- "SCT" # Set default assay to something else than predicted_ADT
adt_data <- breast@assays[["predicted_ADT"]]@data # Extract the data matrix from the predicted_ADT assay
gene_names <- rownames(adt_data) # Get current gene names
new_gene_names <- sapply(gene_names, function(g) adt_map[[g]] %||% g) # Map the gene names using adt_map
rownames(adt_data) <- new_gene_names # Update with protein names
breast[["predicted_ADT"]] <- NULL # Remove the old predicted_ADT assay
breast[["predicted_ADT"]] <- CreateAssayObject(data = adt_data, key = "predictedadt_") # Add the new assay with correct names
saveRDS(breast, "./local/spatial/breast/breast_predicted_protein_names_predicted_ADT.rds")


# Read RDS with protein names
breast <- readRDS("./local/spatial/breast/breast_predicted_protein_names_predicted_ADT.rds")

# Set default assay to be predicted_ADT
DefaultAssay(breast) <- "predicted_ADT"

# Make list of the protein names
proteins <- unique(rownames(breast@assays[["predicted_ADT"]]@data))

# Plot predicted ADT values (spatial) for each protein
for (protein in proteins) {
  SpatialFeaturePlot(breast, features = c(protein))
  ggsave(paste0("./local/spatial/breast/figures/predicted_ADT/", protein, "_breast_spatial.pdf"), width = 20/2, height = 29/2, units = 'cm')
  ggsave(paste0("./local/spatial/breast/figures/predicted_ADT/", protein, "_breast_spatial.png"), width = 20/2, height = 29/2, units = 'cm')
}

# Plot predicted ADt values (spatial) for each protein together with cell annotation
for (protein in proteins) {
  p1 <- SpatialFeaturePlot(breast, features = c(protein))
  p2 <- SpatialDimPlot(breast, group.by = "cell_type", images = "slice1")
  p1 + p2
  ggsave(paste0("./local/spatial/breast/figures/predicted_ADT/annotation_", protein, "_breast_spatial.pdf"), width = 15, height = 10, units = 'cm')
  ggsave(paste0("./local/spatial/breast/figures/predicted_ADT/annotation_", protein, "_breast_spatial.png"), width = 15, height = 10, units = 'cm')
}

