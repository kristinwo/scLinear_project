library(scLinear)
set.seed(42)
options(timeout = 200)

# File: "Gene expression - Feature / cell matrix (filtered)"

# Download the cell matrix file into the local directory and untar it
url <- "https://cf.10xgenomics.com/samples/cell-exp/7.2.0/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Breast_Cancer_BC3_AB3/4plex_DTC_kidney_lung_breast_TotalSeqC_multiplex_Breast_Cancer_BC3_AB3_count_sample_filtered_feature_bc_matrix.tar.gz"
destfile <-"local/breast/Breast_Cancer_BC3_AB3_count_sample_filtered_feature_bc_matrix.tar.gz"   
download.file(url, destfile)
untar(destfile, exdir = "local/breast")

# Create a Seurat object from the data
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

## Create a training and a test set
set.seed(42)
indx <- sample(1:length(colnames(breast)), size = length(colnames(breast)), replace = FALSE)
breast_train <- breast[,indx[1:4000]]
breast_test <- breast[,indx[4001:length(colnames(breast))]]
saveRDS(breast_train ,"./local/breast/breast_train.rds")
saveRDS(breast_test ,"./local/breast/breast_test.rds")

## create predictor
pipe <- create_adt_predictor()

## train predictor
pipe <- fit_predictor(pipe = pipe,
                      gexp_train = breast_train@assays[["RNA"]],
                      adt_train = breast_train@assays[["ADT"]],
                      normalize_gex = TRUE,
                      normalize_adt = TRUE)

## save the trained model
save_trained_model(pipe = pipe, file = "./local/breast/trained_model.joblib")
#> NULL

# load the trained model
pipe <- create_adt_predictor()
pipe <- load_pretrained_model(pipe, file = "./local/breast/trained_model.joblib")

## evaluate predictor
eval_res <- evaluate_predictor(pipe = pipe,
                               gexp_test = breast_test@assays[["RNA"]],
                               adt_test = breast_test@assays[["ADT"]],
                               normalize_gex = TRUE,
                               normalize_adt = TRUE)
#> RMSE: 0.3568380347745682
#> Pearson correlation: 0.8707220850274892
#> Spearman correlation: 0.7528749146936792

print(eval_res)
#>        RMSE   Pearson  Spearman
#> 1 0.356838 0.8707221 0.7528749

## add the predicted adt assay
breast_test[["predicted_ADT"]] <-  adt_predict(pipe = pipe,
                                             gexp = breast_test@assays[["RNA"]],
                                             normalize = TRUE)

saveRDS(breast ,"./local/breast/breast_predicted.rds")


# Train a new model for each cell type ------------------------------------

# Load train and test data
breast_train <- readRDS("./local/breast/breast_train.rds")
breast_test <- readRDS("./local/breast/breast_test.rds")

# Find the common cell types for both datasets
cell_types_train <- unique(breast_train@meta.data[["cell_type"]])
cell_types_test <- unique(breast_test@meta.data[["cell_type"]])
cell_types <- intersect(cell_types_train, cell_types_test)

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
}

#[1] "Duct epithelial cell"
#       RMSE  Pearson  Spearman
#1 0.3569798 0.870793 0.7531325

#[1] "Fibroblast"
#       RMSE   Pearson  Spearman
#1 0.3570126 0.8706443 0.7529327

#[1] "Monocyte"
#       RMSE   Pearson  Spearman
#1 0.3568552 0.8707682 0.7531734

#[1] "Meso-epithelial cell"
#       RMSE   Pearson Spearman
#1 0.3569257 0.8707558 0.753126