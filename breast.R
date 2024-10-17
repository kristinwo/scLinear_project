library(scLinear)
set.seed(42)
options(timeout = 200)

# File: "Gene expression - Feature / cell matrix (filtered)"

# Download the cell matrix file into the local directory and untar it
url <- "https://cf.10xgenomics.com/samples/cell-exp/7.2.0/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Lung_Cancer_BC4_AB4/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Lung_Cancer_BC4_AB4_count_sample_filtered_feature_bc_matrix.tar.gz"
destfile <-"local/lung/Lung_Cancer_BC4_AB4_count_sample_filtered_feature_bc_matrix.tar.gz"   
download.file(url, destfile)
untar(destfile, exdir = "local/lung")

# Create a Seurat object from the data
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

## Create a training and a test set
set.seed(42)
indx <- sample(1:length(colnames(lung)), size = length(colnames(lung)), replace = FALSE)
lung_train <- lung[,indx[1:4000]]
lung_test <- lung[,indx[4001:length(colnames(lung))]]
saveRDS(lung_train ,"./local/lung/lung_train.rds")
saveRDS(lung_test ,"./local/lung/lung_test.rds")

## create predictor
pipe <- create_adt_predictor()

## train predictor
pipe <- fit_predictor(pipe = pipe,
                      gexp_train = lung_train@assays[["RNA"]],
                      adt_train = lung_train@assays[["ADT"]],
                      normalize_gex = TRUE,
                      normalize_adt = TRUE)

## save the trained model
save_trained_model(pipe = pipe, file = "./local/lung/trained_model.joblib")
#> NULL

# load the trained model
pipe <- create_adt_predictor()
pipe <- load_pretrained_model(pipe, file = "./local/lung/trained_model.joblib")

## evaluate predictor
eval_res <- evaluate_predictor(pipe = pipe,
                               gexp_test = lung_test@assays[["RNA"]],
                               adt_test = lung_test@assays[["ADT"]],
                               normalize_gex = TRUE,
                               normalize_adt = TRUE)
#> RMSE: 0.41868590614209306
#> Pearson correlation: 0.9372925884666689
#> Spearman correlation: 0.8452200898462052

print(eval_res)
#>        RMSE   Pearson  Spearman
#> 1 0.4186859 0.9372926 0.8452201

## add the predicted adt assay
lung_test[["predicted_ADT"]] <-  adt_predict(pipe = pipe,
                                             gexp = lung_test@assays[["RNA"]],
                                             normalize = TRUE)

saveRDS(lung ,"./local/lung/lung_predicted.rds")


# Train a new model for each cell type ------------------------------------

# Load train and test data
lung_train <- readRDS("./local/lung/lung_train.rds")
lung_test <- readRDS("./local/lung/lung_test.rds")

# Find the common cell types for both datasets
cell_types_train <- unique(lung_train@meta.data[["cell_type"]])
cell_types_test <- unique(lung_test@meta.data[["cell_type"]])
cell_types <- intersect(cell_types_train, cell_types_test)

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
}

#[1] "Lung epithelial cell"
#       RMSE   Pearson Spearman
#1 0.4189266 0.9372853 0.845312

#[1] "Monocyte"
#       RMSE   Pearson  Spearman
#1 0.4182823 0.9374986 0.8450287

#[1] "T"
#      RMSE   Pearson  Spearman
#1 0.418664 0.9373326 0.8462056

#[1] "Mast cells"
#       RMSE   Pearson  Spearman
#1 0.4190389 0.9371972 0.8453142

#[1] "Meso-epithelial cell"
#       RMSE   Pearson  Spearman
#1 0.4185102 0.9372607 0.8452336

#[1] "B"
#       RMSE   Pearson Spearman
#1 0.4183284 0.9373249 0.845845

#[1] "Fibroblast"
#      RMSE   Pearson  Spearman
#1 0.418867 0.9371014 0.8451156

