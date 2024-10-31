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
url <- "https://cf.10xgenomics.com/samples/cell-exp/7.2.0/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Lung_Cancer_BC4_AB4/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Lung_Cancer_BC4_AB4_count_sample_filtered_feature_bc_matrix.tar.gz"
destfile <-"local/lung/Lung_Cancer_BC4_AB4_count_sample_filtered_feature_bc_matrix.tar.gz"   
download.file(url, destfile)
untar(destfile, exdir = "local/lung")


# Create Seurat object ----------------------------------------------------

data_dir <- "local/lung/sample_filtered_feature_bc_matrix"
lung.data <- Seurat::Read10X(data.dir = data_dir)
rownames(x = lung.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "", x = rownames(x = lung.data[["Antibody Capture"]]))
lung <- Seurat::CreateSeuratObject(counts = lung.data[["Gene Expression"]], min.cells = 1, min.features = 1)
lung[["ADT"]] <- Seurat::CreateAssayObject(lung.data[["Antibody Capture"]][, colnames(x = lung)])
Seurat::DefaultAssay(lung) <- "RNA"

saveRDS(lung, "./local/lung/lung.rds")


# Prepare data ------------------------------------------------------------

lung <- readRDS("./local/lung/lung.rds")
lung <- prepare_data(lung,
                     integrate_data = FALSE,
                     annotation_selfCluster = TRUE, 
                     remove_empty_droplets = FALSE)

saveRDS(lung ,"./local/lung/lung_prepared.rds")


# Train a new model -------------------------------------------------------

lung <- readRDS("./local/lung/lung_prepared.rds")

# Create a training and a test set
set.seed(42)
indx <- sample(1:length(colnames(lung)), size = length(colnames(lung)), replace = FALSE)
lung_train <- lung[,indx[1:4000]]
lung_test <- lung[,indx[4001:length(colnames(lung))]]
saveRDS(lung_train ,"./local/lung/lung_train.rds")
saveRDS(lung_test ,"./local/lung/lung_test.rds")

# Load training and test set
lung_train <- readRDS("./local/lung/lung_train.rds")
lung_test <- readRDS("./local/lung/lung_test.rds")

#Create predictor
pipe <- create_adt_predictor()

# Train predictor
pipe <- fit_predictor(pipe = pipe,
                      gexp_train = lung_train@assays[["RNA"]],
                      adt_train = lung_train@assays[["ADT"]],
                      normalize_gex = TRUE,
                      normalize_adt = TRUE)

# Save the trained model
save_trained_model(pipe = pipe, file = "./local/lung/trained_model.joblib")

# Load the trained model
pipe <- create_adt_predictor()
pipe <- load_pretrained_model(pipe, file = "./local/lung/trained_model.joblib")

# Evaluate predictor
eval_res <- evaluate_predictor(pipe = pipe,
                               gexp_test = lung_test@assays[["RNA"]],
                               adt_test = lung_test@assays[["ADT"]],
                               normalize_gex = TRUE,
                               normalize_adt = TRUE)

print(eval_res)

# Save evaluation results as a csv
write.table(eval_res, "./local/lung/evaluation_results_all_lung.csv")

# Add the predicted adt assay
lung_test[["predicted_ADT"]] <-  adt_predict(pipe = pipe,
                                             gexp = lung_test@assays[["RNA"]],
                                             normalize = TRUE)
# Normalize ADT values
lung_test <- Seurat::NormalizeData(lung_test, normalization.method = "CLR", margin = 2, assay = "ADT")

saveRDS(lung_test ,"./local/lung/lung_predicted.rds")


# Train a new model for each cell type ------------------------------------

# Load train and test data
lung_train <- readRDS("./local/lung/lung_train.rds")
lung_test <- readRDS("./local/lung/lung_test.rds")

# Find the common cell types for both datasets
cell_types_train <- unique(lung_train@meta.data[["cell_type"]])
cell_types_test <- unique(lung_test@meta.data[["cell_type"]])
cell_types <- intersect(cell_types_train, cell_types_test)

# Create a data frame to store the evaluation results for each cell type
results <- data.frame(RMSE = numeric(),
                      Pearson = numeric(),
                      Spearman = numeric(),
                      cell_type = character(),
                      stringsAsFactors = FALSE)


for (cell_type in cell_types) {
  
  # Subset the train and and test data for the specific cell type
  lung_train_subset <- subset(lung_train, subset = lung_train@meta.data[["cell_type"]] == cell_type)
  lung_test_subset <- subset(lung_test, subset = lung_test@meta.data[["cell_type"]] == cell_type)
  
  # Create predictor
  pipe <- create_adt_predictor()
  
  # Train predictor
  pipe <- fit_predictor(pipe = pipe,
                        gexp_train = lung_train_subset@assays[["RNA"]],
                        adt_train = lung_train_subset@assays[["ADT"]],
                        normalize_gex = TRUE,
                        normalize_adt = TRUE)
  
  # Save the trained model
  filename <- paste0("./local/lung/trained_model_", cell_type, ".joblib")
  save_trained_model(pipe = pipe, file = filename)
  
  # Load the trained model
  pipe <- create_adt_predictor()
  pipe <- load_pretrained_model(pipe, file = filename)
  
  # Evaluate predictor
  eval_res <- evaluate_predictor(pipe = pipe,
                                 gexp_test = lung_test_subset@assays[["RNA"]],
                                 adt_test = lung_test_subset@assays[["ADT"]],
                                 normalize_gex = TRUE,
                                 normalize_adt = TRUE)
  
  print(cell_type)
  print(eval_res)
  
  
  
  # Add the predicted adt assay
  lung_test_subset[["predicted_ADT"]] <-  adt_predict(pipe = pipe,
                                                        gexp = lung_test_subset@assays[["RNA"]],
                                                        normalize = TRUE)
  
  # Add the evaluation results to the results data frame
  results <- rbind(results, data.frame(RMSE = eval_res$RMSE,
                                       Pearson = eval_res$Pearson,
                                       Spearman = eval_res$Spearman,
                                       cell_type = cell_type))
}

# Write the evaluation results to a csv file
write.table(results, file = "./local/lung/evaluation_results_celltype_lung.csv", 
            sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Create correlation figure for each protein ------------------------------

lung <- readRDS("./local/lung/lung_predicted.rds")

# Create data frame of real and predicted ADT values
adt_real <- as.matrix(lung@assays[["ADT"]]@data)
adt_real <- adt_real %>% as.data.frame() %>% rownames_to_column("gene") %>% pivot_longer(!gene,names_to = "cell", values_to = "real")
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
adt_real$gene <- adt_map[adt_real$gene]
adt_predicted_sclinear$gene <- adt_map[adt_predicted_sclinear$gene]

# Get cell types for each cell
meta <- lung@meta.data %>% rownames_to_column("cell") %>% dplyr::select("cell", "cell_type")

# Create data frame of real and predicted ADT values
DF <- adt_real %>% full_join(adt_predicted_sclinear, by = c("gene", "cell")) %>% full_join(meta, by = c("cell"))
DF <- DF %>% arrange(gene)
DF$gene <- factor(DF$gene, levels = unique(DF$gene))

write.table(DF, file = "./local/lung/bit_table_lung.csv", sep = ",", col.names = TRUE, row.names = FALSE)

lung <- read.table("./local/lung/bit_table_lung.csv", header = T, sep=',')
lung <- lung %>% distinct()

# Calculate pearson correlation for each protein
cors_pearson <- lung %>% group_by(gene) %>% dplyr::summarise(correlation = cor(real, predicted_sclinear, method="pearson"))
write.csv(cors_pearson, "./local/lung/pearson_correlation_proteins_lung.csv")

# Calculate spearman correlation for each protein
cors_spearman <- lung %>% group_by(gene) %>% dplyr::summarise(correlation = cor(real, predicted_sclinear, method="spearman"))
write.csv(cors_spearman, "./local/lung/spearman_correlation_proteins_lung.csv")

# Calculate RMSE for each protein
rmse <- lung %>% group_by(gene) %>% dplyr::summarise(rmse = rmse(real, predicted_sclinear))
write.csv(cors_spearman, "./local/lung/rmse_proteins_lung.csv")


# Plot the pearson correlation for each protein
ggplot(cors_pearson, aes(x=reorder(gene, correlation), y= correlation,fill=correlation)) + 
  geom_bar(stat="identity", col="black") + 
  coord_flip() +
  theme_classic() + 
  scale_fill_gradientn(colours = inferno(11)) + 
  ylab("Pearson\n(Real vs ScLinear)") + 
  xlab("Protein") + 
  theme(legend.position = "none") +
  ggtitle("Lung CITE-seq") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y = 1)

ggsave("./local/lung/figures/lung_Pearson_proteins.pdf", width = 20/2, height = 29/2, units = 'cm')
ggsave("./local/lung/figures/lung_Pearson_proteins.png", width = 20/2, height = 29/2, units = 'cm')

# Plot the spearman correlation for each protein
ggplot(cors_spearman, aes(x=reorder(gene, correlation), y= correlation,fill=correlation)) + 
  geom_bar(stat="identity", col="black") + 
  coord_flip() +
  theme_classic() + 
  scale_fill_gradientn(colours = inferno(11)) + 
  ylab("Spearman\n(Real vs ScLinear)") + 
  xlab("Protein") + 
  theme(legend.position = "none") +
  ggtitle("Lung CITE-seq") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y = 1)

ggsave("./local/lung/figures/lung_Spearman_proteins.pdf", width = 20/2, height = 29/2, units = 'cm')
ggsave("./local/lung/figures/lung_Spearman_proteins.png", width = 20/2, height = 29/2, units = 'cm')

# Plot the RMSE for each protein
ggplot(rmse, aes(x=reorder(gene, rmse), y= rmse,fill=rmse)) + 
  geom_bar(stat="identity", col="black") + 
  coord_flip() +
  theme_classic() + 
  scale_fill_gradientn(colours = inferno(11)) + 
  ylab("RMSE\n(Real vs ScLinear)") + 
  xlab("Protein") + 
  theme(legend.position = "none") +
  ggtitle("Lung CITE-seq") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y = 1)

ggsave("./local/lung/figures/lung_RMSE_proteins.pdf", width = 20/2, height = 29/2, units = 'cm')
ggsave("./local/lung/figures/lung_RMSE_proteins.png", width = 20/2, height = 29/2, units = 'cm')
