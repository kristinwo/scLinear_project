library(scLinear)
library(tidyverse)
library(ggplot2)
library(Seurat)
library(svglite)
library(Metrics)
library(pals)
library(ggpubr)
library(Matrix)
library(grid)
library(ggrepel)
library(hdf5r)
source("helper_functions.R")
set.seed(42)
options(timeout = 200)


# Get data ----------------------------------------------------------------

# Download the cell matrix file into the local directory and untar it

# Rep 1
url1 <- "https://cf.10xgenomics.com/samples/cell-exp/7.2.0/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Kidney_Cancer1_BC1_AB1/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Kidney_Cancer1_BC1_AB1_count_sample_filtered_feature_bc_matrix.tar.gz"
destfile1 <-"local/kidney/kidney_rep1_filtered_feature_bc_matrix.tar.gz"   
download.file(url1, destfile1)
untar(destfile1, exdir = "local/kidney")
file.rename("./local/kidney/sample_filtered_feature_bc_matrix", "./local/kidney/kidney_rep1_sample_filtered_feature_bc_matrix")


# Rep 2
url2 <- "https://cf.10xgenomics.com/samples/cell-exp/7.2.0/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Kidney_Cancer2_BC2_AB2/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Kidney_Cancer2_BC2_AB2_count_sample_filtered_feature_bc_matrix.tar.gz"
destfile2 <-"local/kidney/kidney_rep2_filtered_feature_bc_matrix.tar.gz"   
download.file(url2, destfile2)
untar(destfile2, exdir = "local/kidney")
file.rename("local/kidney/sample_filtered_feature_bc_matrix", "local/kidney/kidney_rep2_sample_filtered_feature_bc_matrix")


# Create Seurat object ----------------------------------------------------

# Rep 1
data_dir <- "local/kidney/kidney_rep1_sample_filtered_feature_bc_matrix"
kidneyRep1.data <- Seurat::Read10X(data.dir = data_dir)
rownames(x = kidneyRep1.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "", x = rownames(x = kidneyRep1.data[["Antibody Capture"]]))
kidneyRep1 <- Seurat::CreateSeuratObject(counts = kidneyRep1.data[["Gene Expression"]], min.cells = 1, min.features = 1)
kidneyRep1[["ADT"]] <- Seurat::CreateAssayObject(kidneyRep1.data[["Antibody Capture"]][, colnames(x = kidneyRep1)])
Seurat::DefaultAssay(kidneyRep1) <- "RNA"

saveRDS(kidneyRep1, "./local/kidney/kidneyRep1.rds")

# Rep 2
data_dir <- "local/kidney/kidney_rep2_sample_filtered_feature_bc_matrix"
kidneyRep2.data <- Seurat::Read10X(data.dir = data_dir)
rownames(x = kidneyRep2.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "", x = rownames(x = kidneyRep2.data[["Antibody Capture"]]))
kidneyRep2 <- Seurat::CreateSeuratObject(counts = kidneyRep2.data[["Gene Expression"]], min.cells = 1, min.features = 1)
kidneyRep2[["ADT"]] <- Seurat::CreateAssayObject(kidneyRep2.data[["Antibody Capture"]][, colnames(x = kidneyRep2)])
Seurat::DefaultAssay(kidneyRep2) <- "RNA"

saveRDS(kidneyRep2, "./local/kidney/kidneyRep2.rds")


# Prepare data ------------------------------------------------------------

# Rep 1
kidneyRep1 <- readRDS("./local/kidney/kidneyRep1.rds")
kidneyRep1 <- prepare_data(kidneyRep1,
                           integrate_data = FALSE,
                           annotation_selfCluster = TRUE, 
                           remove_empty_droplets = FALSE)


saveRDS(kidneyRep1 ,"./local/kidney/kidneyRep1_prepared.rds")

# Rep 2
kidneyRep2 <- readRDS("./local/kidney/kidneyRep2.rds")
kidneyRep2 <- prepare_data(kidneyRep2,
                           integrate_data = FALSE,
                           annotation_selfCluster = TRUE, 
                           remove_empty_droplets = FALSE)


saveRDS(kidneyRep2 ,"./local/kidney/kidneyRep2_prepared.rds")



# Train a new model -------------------------------------------------------

# Create a training and a test set
kidney_train <- readRDS("./local/kidney/kidneyRep1_prepared.rds")
kidney_test <- readRDS("./local/kidney/kidneyRep2_prepared.rds")

# Create predictor
pipe <- create_adt_predictor()

# Train predictor
pipe <- fit_predictor(pipe = pipe,
                      gexp_train = kidney_train@assays[["RNA"]],
                      adt_train = kidney_train@assays[["ADT"]],
                      normalize_gex = TRUE,
                      normalize_adt = TRUE)

# Save the trained model
save_trained_model(pipe = pipe, file = "./local/kidney/trained_model.joblib")
#> NULL

# Load the trained model
pipe <- create_adt_predictor()
pipe <- load_pretrained_model(pipe, file = "./local/kidney/trained_model.joblib")

# Evaluate predictor
eval_res <- evaluate_predictor(pipe = pipe,
                               gexp_test = kidney_test@assays[["RNA"]],
                               adt_test = kidney_test@assays[["ADT"]],
                               normalize_gex = TRUE,
                               normalize_adt = TRUE)

print(eval_res)

# Save evaluation results as a csv
write.table(eval_res, "./local/kidney/evaluation_results_all_kidney.csv")

# Add the predicted adt assay
kidney_test[["predicted_ADT"]] <-  adt_predict(pipe = pipe,
                                               gexp = kidney_test@assays[["RNA"]],
                                               normalize = TRUE)

# Normalize ADT values
kidney_test <- Seurat::NormalizeData(kidney_test, normalization.method = "CLR", margin = 2, assay = "ADT")

saveRDS(kidney_test,"./local/kidney/kidney_predicted.rds")

# Train a new model for each cell type ------------------------------------

# Load train and test data
kidney_train <- readRDS("./local/kidney/kidneyRep1_prepared.rds")
kidney_test <- readRDS("./local/kidney/kidneyRep2_prepared.rds")

# Find the common cell types for both datasets
cell_types_train <- unique(kidney_train@meta.data[["cell_type"]])
cell_types_test <- unique(kidney_test@meta.data[["cell_type"]])
cell_types <- intersect(cell_types_train, cell_types_test)

# Create a data frame to store the evaluation results for each cell type
results <- data.frame(RMSE = numeric(),
                      Pearson = numeric(),
                      Spearman = numeric(),
                      cell_type = character(),
                      stringsAsFactors = FALSE)


for (cell_type in cell_types) {
  
  # Subset the train and and test data for the specific cell type
  kidney_train_subset <- subset(kidney_train, subset = kidney_train@meta.data[["cell_type"]] == cell_type)
  kidney_test_subset <- subset(kidney_test, subset = kidney_test@meta.data[["cell_type"]] == cell_type)
  
  # Create predictor
  pipe <- create_adt_predictor()
  
  # Train predictor
  pipe <- fit_predictor(pipe = pipe,
                        gexp_train = kidney_train_subset@assays[["RNA"]],
                        adt_train = kidney_train_subset@assays[["ADT"]],
                        normalize_gex = TRUE,
                        normalize_adt = TRUE)
  
  # Save the trained model
  filename <- paste0("./local/kidney/trained_model_", cell_type, ".joblib")
  save_trained_model(pipe = pipe, file = filename)
  
  # Load the trained model
  pipe <- create_adt_predictor()
  pipe <- load_pretrained_model(pipe, file = filename)
  
  # Evaluate predictor
  eval_res <- evaluate_predictor(pipe = pipe,
                                 gexp_test = kidney_test_subset@assays[["RNA"]],
                                 adt_test = kidney_test_subset@assays[["ADT"]],
                                 normalize_gex = TRUE,
                                 normalize_adt = TRUE)
  
  print(cell_type)
  print(eval_res)
  
  
  
  # Add the predicted adt assay
  kidney_test_subset[["predicted_ADT"]] <-  adt_predict(pipe = pipe,
                                                        gexp = kidney_test_subset@assays[["RNA"]],
                                                        normalize = TRUE)
  
  # Add the evaluation results to the results data frame
  results <- rbind(results, data.frame(RMSE = eval_res$RMSE,
                                       Pearson = eval_res$Pearson,
                                       Spearman = eval_res$Spearman,
                                       cell_type = cell_type))
}

# Write the evaluation results to a csv file
write.table(results, file = "./local/kidney/evaluation_results_celltype_kidney.csv", 
            sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Create correlation figure for each protein ------------------------------

kidney <- readRDS("./local/kidney/kidney_predicted.rds")

# Create data frame of real and predicted ADT values
adt_real <- as.matrix(kidney@assays[["ADT"]]@data)
adt_real <- adt_real %>% as.data.frame() %>% rownames_to_column("gene") %>% pivot_longer(!gene,names_to = "cell", values_to = "real")
adt_predicted_sclinear <- as.matrix(kidney@assays[["predicted_ADT"]]@data)
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
adt_real$gene <- adt_map[adt_real$gene]
adt_predicted_sclinear$gene <- adt_map[adt_predicted_sclinear$gene]

# Get cell types for each cell
meta <- kidney@meta.data %>% rownames_to_column("cell") %>% dplyr::select("cell", "cell_type")

# Create data frame of real and predicted ADT values
DF <- adt_real %>% full_join(adt_predicted_sclinear, by = c("gene", "cell")) %>% full_join(meta, by = c("cell"))
DF <- DF %>% arrange(gene)
DF$gene <- factor(DF$gene, levels = unique(DF$gene))

# Write bit table
write.table(DF, file = "./local/kidney/bit_table_kidney.csv", sep = ",", col.names = TRUE, row.names = FALSE)

# Read bit table
kidney <- read.table("./local/kidney/bit_table_kidney.csv", header = T, sep=',')
kidney <- kidney %>% distinct()

# Calculate pearson correlation for each protein
cors_pearson <- kidney %>% group_by(gene) %>% dplyr::summarise(correlation = cor(real, predicted_sclinear, method="pearson"))
write.csv(cors_pearson, "./local/kidney/pearson_correlation_proteins_kidney.csv")

# Calculate spearman correlation for each protein
cors_spearman <- kidney %>% group_by(gene) %>% dplyr::summarise(correlation = cor(real, predicted_sclinear, method="spearman"))
write.csv(cors_spearman, "./local/kidney/spearman_correlation_proteins_kidney.csv")

# Calculate RMSE for each protein
rmse <- kidney %>% group_by(gene) %>% dplyr::summarise(rmse = rmse(real, predicted_sclinear))
write.csv(cors_spearman, "./local/kidney/rmse_proteins_kidney.csv")


# Plot the pearson correlation for each protein
ggplot(cors_pearson, aes(x=reorder(gene, correlation), y= correlation,fill=correlation)) + 
  geom_bar(stat="identity", col="black") + 
  coord_flip() +
  theme_classic() + 
  scale_fill_gradientn(colours = inferno(11)) + 
  ylab("Pearson\n(Real vs ScLinear)") + 
  xlab("Protein") + 
  theme(legend.position = "none") +
  ggtitle("Kidney CITE-seq") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y = 1)

ggsave("./local/kidney/figures/Kidney_Pearson_proteins.pdf", width = 20/2, height = 29/2, units = 'cm')
ggsave("./local/kidney/figures/Kidney_Pearson_proteins.png", width = 20/2, height = 29/2, units = 'cm')

# Plot the spearman correlation for each protein
ggplot(cors_spearman, aes(x=reorder(gene, correlation), y= correlation,fill=correlation)) + 
  geom_bar(stat="identity", col="black") + 
  coord_flip() +
  theme_classic() + 
  scale_fill_gradientn(colours = inferno(11)) + 
  ylab("Spearman\n(Real vs ScLinear)") + 
  xlab("Protein") + 
  theme(legend.position = "none") +
  ggtitle("Kidney CITE-seq") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y = 1)

ggsave("./local/kidney/figures/Kidney_Spearman_proteins.pdf", width = 20/2, height = 29/2, units = 'cm')
ggsave("./local/kidney/figures/Kidney_Spearman_proteins.png", width = 20/2, height = 29/2, units = 'cm')

# Plot the RMSE for each protein
ggplot(rmse, aes(x=reorder(gene, rmse), y= rmse,fill=rmse)) + 
  geom_bar(stat="identity", col="black") + 
  coord_flip() +
  theme_classic() + 
  scale_fill_gradientn(colours = inferno(11)) + 
  ylab("RMSE\n(Real vs ScLinear)") + 
  xlab("Protein") + 
  theme(legend.position = "none") +
  ggtitle("Kidney CITE-seq") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y = 1)

ggsave("./local/kidney/figures/Kidney_RMSE_proteins.pdf", width = 20/2, height = 29/2, units = 'cm')
ggsave("./local/kidney/figures/Kidney_RMSE_proteins.png", width = 20/2, height = 29/2, units = 'cm')


