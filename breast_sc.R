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
url <- "https://cf.10xgenomics.com/samples/cell-exp/7.2.0/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Breast_Cancer_BC3_AB3/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Breast_Cancer_BC3_AB3_count_sample_filtered_feature_bc_matrix.tar.gz"
destfile <-"local/breast/Breast_Cancer_BC3_AB3_count_sample_filtered_feature_bc_matrix.tar.gz"   
download.file(url, destfile)
untar(destfile, exdir = "local/breast")


# Create Seurat object ----------------------------------------------------

data_dir <- "local/breast/sample_filtered_feature_bc_matrix"
breast.data <- Seurat::Read10X(data.dir = data_dir)
rownames(x = breast.data[["Antibody Capture"]]) <- gsub(pattern = "_[control_]*TotalSeqB", replacement = "", x = rownames(x = breast.data[["Antibody Capture"]]))
breast <- Seurat::CreateSeuratObject(counts = breast.data[["Gene Expression"]], min.cells = 1, min.features = 1)
breast[["ADT"]] <- Seurat::CreateAssayObject(breast.data[["Antibody Capture"]][, colnames(x = breast)])
Seurat::DefaultAssay(breast) <- "RNA"

saveRDS(breast, "./local/breast/breast.rds")


# Prepare data ------------------------------------------------------------

breast <- readRDS("./local/breast/breast.rds")
breast <- prepare_data(breast,
                     integrate_data = FALSE,
                     annotation_selfCluster = TRUE, 
                     remove_empty_droplets = FALSE)

saveRDS(breast ,"./local/breast/breast_prepared.rds")


# Train a new model -------------------------------------------------------

breast <- readRDS("./local/breast/breast_prepared.rds")

# Create a training and a test set
set.seed(42)
indx <- sample(1:length(colnames(breast)), size = length(colnames(breast)), replace = FALSE)
breast_train <- breast[,indx[1:4000]]
breast_test <- breast[,indx[4001:length(colnames(breast))]]
saveRDS(breast_train ,"./local/breast/breast_train.rds")
saveRDS(breast_test ,"./local/breast/breast_test.rds")

# Load training and test set
breast_train <- readRDS("./local/breast/breast_train.rds")
breast_test <- readRDS("./local/breast/breast_test.rds")

# Create predictor
pipe <- create_adt_predictor()

# Train predictor
pipe <- fit_predictor(pipe = pipe,
                      gexp_train = breast_train@assays[["RNA"]],
                      adt_train = breast_train@assays[["ADT"]],
                      normalize_gex = TRUE,
                      normalize_adt = TRUE)

# Save the trained model
save_trained_model(pipe = pipe, file = "./local/breast/trained_model.joblib")

# Load the trained model
pipe <- create_adt_predictor()
pipe <- load_pretrained_model(pipe, file = "./local/breast/trained_model.joblib")

# Evaluate predictor
eval_res <- evaluate_predictor(pipe = pipe,
                               gexp_test = breast_test@assays[["RNA"]],
                               adt_test = breast_test@assays[["ADT"]],
                               normalize_gex = TRUE,
                               normalize_adt = TRUE)

print(eval_res)

# Save evaluation results as a csv
write.table(eval_res, "./local/breast/evaluation_results_all_breast.csv")

# Add the predicted adt assay
breast_test[["predicted_ADT"]] <- adt_predict(pipe = pipe,
                                              gexp = breast_test@assays[["RNA"]],
                                              normalize = TRUE)

# Normalize ADT values
breast_test <- Seurat::NormalizeData(breast_test, normalization.method = "CLR", margin = 2, assay = "ADT")


saveRDS(breast_test,"./local/breast/breast_predicted.rds")


# Train a new model for each cell type ------------------------------------

# Load train and test data
breast_train <- readRDS("./local/breast/breast_train.rds")
breast_test <- readRDS("./local/breast/breast_test.rds")

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
  filename <- paste0("./local/breast/trained_model_", cell_type, ".joblib")
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
write.table(results, file = "./local/breast/evaluation_results_celltype_breast.csv", 
            sep = ",", row.names = FALSE, col.names = TRUE, quote = FALSE)


# Create correlation figure for each protein ------------------------------

breast <- readRDS("./local/breast/breast_predicted.rds")

# Create data frame of real and predicted ADT values
adt_real <- as.matrix(breast@assays[["ADT"]]@data)
adt_real <- adt_real %>% as.data.frame() %>% rownames_to_column("gene") %>% pivot_longer(!gene,names_to = "cell", values_to = "real")
adt_predicted_sclinear <- as.matrix(breast@assays[["predicted_ADT"]]@data)
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
meta <- breast@meta.data %>% rownames_to_column("cell") %>% dplyr::select("cell", "cell_type")

# Create data frame of real and predicted ADT values
DF <- adt_real %>% full_join(adt_predicted_sclinear, by = c("gene", "cell")) %>% full_join(meta, by = c("cell"))
DF <- DF %>% arrange(gene)
DF$gene <- factor(DF$gene, levels = unique(DF$gene))

write.table(DF, file = "./local/breast/bit_table_breast.csv", sep = ",", col.names = TRUE, row.names = FALSE)

breast <- read.table("./local/breast/bit_table_breast.csv", header = T, sep=',')
breast <- breast %>% distinct()

# Calculate pearson correlation for each protein
cors_pearson <- breast %>% group_by(gene) %>% dplyr::summarise(correlation = cor(real, predicted_sclinear, method="pearson"))
write.csv(cors_pearson, "./local/breast/pearson_correlation_proteins_breast.csv")

# Calculate spearman correlation for each protein
cors_spearman <- breast %>% group_by(gene) %>% dplyr::summarise(correlation = cor(real, predicted_sclinear, method="spearman"))
write.csv(cors_spearman, "./local/breast/spearman_correlation_proteins_breast.csv")

# Calculate RMSE for each protein
rmse <- breast %>% group_by(gene) %>% dplyr::summarise(rmse = rmse(real, predicted_sclinear))
write.csv(cors_spearman, "./local/breast/rmse_proteins_breast.csv")


# Plot the pearson correlation for each protein
ggplot(cors_pearson, aes(x=reorder(gene, correlation), y= correlation,fill=correlation)) + 
  geom_bar(stat="identity", col="black") + 
  coord_flip() +
  theme_classic() + 
  scale_fill_gradientn(colours = inferno(11)) + 
  ylab("Pearson\n(Real vs ScLinear)") + 
  xlab("Protein") + 
  theme(legend.position = "none") +
  ggtitle("Breast CITE-seq") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y = 1)

ggsave("./local/breast/figures/breast_Pearson_proteins.pdf", width = 20/2, height = 29/2, units = 'cm')
ggsave("./local/breast/figures/breast_Pearson_proteins.png", width = 20/2, height = 29/2, units = 'cm')

# Plot the spearman correlation for each protein
ggplot(cors_spearman, aes(x=reorder(gene, correlation), y= correlation,fill=correlation)) + 
  geom_bar(stat="identity", col="black") + 
  coord_flip() +
  theme_classic() + 
  scale_fill_gradientn(colours = inferno(11)) + 
  ylab("Spearman\n(Real vs ScLinear)") + 
  xlab("Protein") + 
  theme(legend.position = "none") +
  ggtitle("Breast CITE-seq") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y = 1)

ggsave("./local/breast/figures/breast_Spearman_proteins.pdf", width = 20/2, height = 29/2, units = 'cm')
ggsave("./local/breast/figures/breast_Spearman_proteins.png", width = 20/2, height = 29/2, units = 'cm')

# Plot the RMSE for each protein
ggplot(rmse, aes(x=reorder(gene, rmse), y= rmse,fill=rmse)) + 
  geom_bar(stat="identity", col="black") + 
  coord_flip() +
  theme_classic() + 
  scale_fill_gradientn(colours = inferno(11)) + 
  ylab("RMSE\n(Real vs ScLinear)") + 
  xlab("Protein") + 
  theme(legend.position = "none") +
  ggtitle("Breast CITE-seq") + 
  theme(plot.title = element_text(hjust = 0.5)) +
  expand_limits(y = 1)

ggsave("./local/breast/figures/breast_RMSE_proteins.pdf", width = 20/2, height = 29/2, units = 'cm')
ggsave("./local/breast/figures/breast_RMSE_proteins.png", width = 20/2, height = 29/2, units = 'cm')
