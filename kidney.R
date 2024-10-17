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
#Warning messages:
#1: Layer ‘data’ is empty 
#2: Layer ‘scale.data’ is empty 
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
#> RMSE: 0.35495989559489205
#> Pearson correlation: 0.933335317306118
#> Spearman correlation: 0.8703162049287475

print(eval_res)
#>      RMSE   Pearson  Spearman
#1 0.3549599 0.9333353 0.8703162

# Add the predicted adt assay
kidney_test[["predicted_ADT"]] <-  adt_predict(pipe = pipe,
                                               gexp = kidney_test@assays[["RNA"]],
                                               normalize = TRUE)

saveRDS(kidney_test,"./local/kidney/kidney_predicted.rds")

# Train a new model for each cell type ------------------------------------

# Load train and test data
kidney_train <- readRDS("./local/kidney/kidneyRep1_prepared.rds")
kidney_test <- readRDS("./local/kidney/kidneyRep2_prepared.rds")

# Find the common cell types for both datasets
cell_types_train <- unique(kidney_train@meta.data[["cell_type"]])
cell_types_test <- unique(kidney_test@meta.data[["cell_type"]])
cell_types <- intersect(cell_types_train, cell_types_test)

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
}

#[1] "Kidney epithelial cell"
#       RMSE   Pearson  Spearman
#1 0.3552711 0.9333183 0.8699878

#[1] "Duct epithelial cell"
#       RMSE   Pearson  Spearman
#1 0.3554204 0.9332646 0.8699997

#[1] "Meso-epithelial cell"
#       RMSE   Pearson  Spearman
#1 0.3552383 0.9333938 0.8703636

#[1] "T"
#       RMSE   Pearson  Spearman
#1 0.3552776 0.9332597 0.8698043

#[1] "Pericytes"
#       RMSE   Pearson  Spearman
#1 0.3551415 0.9333342 0.8703374

#[1] "Monocyte"
#       RMSE   Pearson  Spearman
#1 0.3554118 0.9331661 0.8701949

#[1] "B"
#       RMSE   Pearson  Spearman
#1 0.3553983 0.9333131 0.8700432


# Create correlation figure for each protein ------------------------------

kidney <- readRDS("./local/kidney/kidney_predicted.rds")

adt_real <- as.matrix(kidney@assays[["ADT"]]@data)
adt_real <- adt_real %>% as.data.frame() %>% rownames_to_column("gene") %>% pivot_longer(!gene,names_to = "cell", values_to = "real")
adt_predicted_sclinear <- as.matrix(kidney@assays[["predicted_ADT"]]@data)
adt_predicted_sclinear <- adt_predicted_sclinear %>% as.data.frame() %>% rownames_to_column("gene") %>% pivot_longer(!gene,names_to = "cell", values_to = "predicted_sclinear")

#### raw RNA counts
rna_raw <- as.matrix(GetAssayData(kidney, assay = "RNA", layer = "counts"))
rna_raw <- rna_raw %>% as.data.frame() %>% rownames_to_column("gene") %>% pivot_longer(!gene,names_to = "cell", values_to = "rna_raw")
#### post processing of RNA data to make it comparable to ADT data

#### add normalized RNA (normalized the same way as used for prediction)
rna_normalized <- as.matrix(scLinear::gexp_normalize(GetAssayData(kidney, assay = "RNA", layer = "counts")))
rna_normalized <- rna_normalized %>% as.data.frame() %>% rownames_to_column("gene") %>% pivot_longer(!gene,names_to = "cell", values_to = "rna_normalized")
#### post processing of RNA data to make it comparable to ADT data

adt_real$adt_gene_names_real <- adt_real$gene

adt_predicted_sclinear$adt_gene_names_predicted <- adt_predicted_sclinear$gene

rna_raw$rna_gene_names_raw <- rna_raw$gene

rna_normalized$rna_gene_names_normalized <- rna_normalized$gene




map <- list(
  
  list(gexp = c("FUT4"), adt = c("CD15")),
  list(gexp = c("FCGR3A"), adt = c("CD16")),
  list(gexp = c("NCAM1"), adt = c("CD56")),
  list(gexp = c("IL7R"), adt = c("CD127")),
  list(gexp = c("PDCD1"), adt = c("CD279")),
  list(gexp = c("CD3E", "CD3D", "CD3G"), adt = c("CD3")),
  list(gexp = c("CD8A", "CD8B"), adt = c("CD8")),
  list(gexp = c("PTPRC"), adt = c("CD45RA", "CD45RO"))
)


for (m in map){
  print(paste0("gexp: ", paste0(m$gexp, collapse = ","), "       adt: ", paste0(m$adt, collapse = ",")))
  
  gexp_names <- m$gexp
  adt_names <- m$adt
  
  if((length(gexp_names) == 1) & (length(adt_names) == 1)){
    ## change rna name to adt names
    rna_raw$gene[rna_raw$gene == c(gexp_names)] <- adt_names
    rna_normalized$gene[rna_normalized$gene == c(gexp_names)] <- adt_names
  }else{
    if((length(gexp_names) > 1) & (length(adt_names) == 1)){
      ## map adt name to many gexp names. each gene compared to adt name.
      genes <- gexp_names
      tmp <- adt_real[adt_real$gene == adt_names,]
      adt_real <- adt_real[!(adt_real$gene == adt_names),]
      for (g in genes){
        tmp_2 <- tmp
        tmp_2 $gene <- g
        adt_real <- rbind(adt_real, tmp_2)
      }
      genes <- gexp_names
      tmp <- adt_predicted_sclinear[adt_predicted_sclinear$gene == adt_names,]
      adt_predicted_sclinear <- adt_predicted_sclinear[!(adt_predicted_sclinear$gene == adt_names),]
      for (g in genes){
        tmp_2 <- tmp
        tmp_2 $gene <- g
        adt_predicted_sclinear <- rbind(adt_predicted_sclinear, tmp_2)
      }
      
    }else{
      if((length(gexp_names) == 1) & (length(adt_names) > 1)){
        genes <- adt_names
        tmp <- rna_raw[rna_raw$gene == gexp_names,]
        rna_raw <- rna_raw[!(rna_raw$gene == gexp_names),]
        for (g in genes){
          tmp_2 <- tmp
          tmp_2 $gene <- g
          rna_raw <- rbind(rna_raw, tmp_2)
        }
        
        genes <- adt_names
        tmp <- rna_normalized[rna_normalized$gene == gexp_names,]
        rna_normalized <- rna_normalized[!(rna_normalized$gene == gexp_names),]
        for (g in genes){
          tmp_2 <- tmp
          tmp_2 $gene <- g
          rna_normalized <- rbind(rna_normalized, tmp_2)
        }
        
      }
    }
  }
  
}


## remove not usefull genes
genes_to_keep <- unique(c(unlist(map, recursive = TRUE), adt_real$gene, adt_predicted_sclinear$gene))
rna_raw <- rna_raw %>% dplyr:: filter(gene %in% genes_to_keep)
rna_normalized <- rna_normalized %>% dplyr:: filter(gene %in% genes_to_keep)

meta <- kidney@meta.data %>% rownames_to_column("cell") %>% dplyr::select("cell", "cell_type")

DF <- adt_real %>% full_join(adt_predicted_sclinear, by = c("gene", "cell")) %>%  full_join(rna_normalized, by = c("gene", "cell")) %>% full_join(rna_raw, by = c("gene", "cell")) %>% full_join(meta, by = c("cell"))

DF <- DF %>% arrange(gene)

DF$gene <- factor(DF$gene, levels = unique(DF$gene))

write.table(DF, file = "./local/kidney/bit_table_kidney.csv", sep = ",", col.names = TRUE, row.names = FALSE)


### ADT predicted correlation comparison

kidney <- read.table("./local/kidney/bit_table_kidney.csv", header = T, sep=',')
kidney <- kidney %>% dplyr::select(-c(rna_normalized, rna_raw))
kidney1 <- na.omit(kidney)

# deduplicate the separation of complexes that was used for RNA/ADT comparison
map <- list(  
  list(gexp = c("CD3E", "CD3D", "CD3G"), adt = c("CD3")),
  list(gexp = c("CD8A", "CD8B"), adt = c("CD8"))
)
for (m in map){
  print(paste0("gexp: ", paste0(m$gexp, collapse = ","), "       adt: ", paste0(m$adt, collapse = ",")))
  gexp_names <- m$gexp
  adt_name <- m$adt
  for(gexp_name in gexp_names){
    kidney1$gene[kidney1$gene == gexp_name] <- adt_name
  }
}
kidney1 <- kidney1 %>% distinct()


cors <- kidney1 %>% group_by(gene) %>% dplyr::summarise(correlation = cor(real, predicted_sclinear))
ggplot(cors, aes(x=reorder(gene, correlation), y= correlation,fill=correlation)) + geom_bar(stat="identity", col="black") + coord_flip() +
  theme_classic() + scale_fill_gradientn(colours = inferno(11)) + ylab("Pearson\n(Real vs ScLinear)") + xlab("Protein") + theme(legend.position = "none") +
  ggtitle("kidney1 CITE-seq\n(10k cells)") + theme(plot.title = element_text(hjust = 0.5))

ggsave("./local/kidney/figures/E.kidney_Pearson_omit_na_proteins.pdf", width = 20/3, height = 29/3, units = 'cm')


my_pal <- c("#E7E7E0","#A2A187", "#525240", "#03878F", "#F2B531","#ED5958","#68C6A4", "#113245")

kidney12 <- filter(kidney1, gene %in% c("CD19","CD3","CD14","CD56"))
ggplot(kidney12, aes(x=predicted_sclinear, y=real, col=cell_type_2)) + geom_point(alpha = 0.8, size=0.5) + facet_wrap(~gene, scales = "free") +
  scale_color_manual(values = my_pal[4:8]) + theme_classic2() + theme(legend.title = element_blank()) +
  geom_smooth(method=lm, se=FALSE, col=my_pal[2], alpha = 0.5) + geom_rug() +
  xlab("ScLinear - Predicted") + ylab("Real")

ggsave("./local/kidney1/figures/F.kidney1_Markers.pdf", width = 20/3*2, height = 29/3, units = "cm")