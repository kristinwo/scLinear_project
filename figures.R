library(ggplot2)
library(reshape)
library(dplyr)
library(tidyr)
library(scLinear)
set.seed(42)


# Summary of RMSE, Pearson and Spearman for kidney, lung and breast --------

# Read csv files containing RMSE, Pearson and Spearman correlations
results_kidney <- read.table("./local/kidney/evaluation_results_all_kidney.csv")
results_lung <- read.table("./local/lung/evaluation_results_all_lung.csv")
results_breast <- read.table("./local/breast/evaluation_results_all_breast.csv")

# Add a column specifying the tissue the data is from
results_kidney$Tissue <- "Kidney"
results_lung$Tissue <- "Lung"
results_breast$Tissue <- "Breast"

# Merge into one single data frame
results <- rbind(results_kidney, results_lung, results_breast)

# Reshape the data to long format for easier plotting
long_results <- melt(results, id.vars = "Tissue", 
                     measure.vars = c("RMSE", "Pearson", "Spearman"))

# Define colors for each tissue
colors <- c(
  "Lung" = "#7D7D61", 
  "Kidney" = "#7B90AA",
  "Breast" = "#BF878A"
)


# Create the RMSE plot
rmse_plot <- ggplot(subset(long_results, variable == "RMSE"), 
                    aes(x = Tissue, y = value, fill = Tissue)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.4) +
  labs(title = "RMSE",
       y = "RMSE") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 1)) + 
  scale_fill_manual(name = "Tissue", values = colors)
  

ggsave("./figures/RMSE_all.pdf", width = 10, height = 8, units = 'cm')
ggsave("./figures/RMSE_all.png", width = 10, height = 8, units = 'cm')


# Create the Pearson plot
pearson_plot <- ggplot(subset(long_results, variable == "Pearson"), 
                       aes(x = Tissue, y = value, fill = Tissue)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.4) +
  labs(title = "Pearson",
       y = "Pearson Correlation") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(name = "Tissue", values = colors)

ggsave("./figures/Pearson_all.pdf", width = 10, height = 8, units = 'cm')
ggsave("./figures/Pearson_all.png", width = 10, height = 8, units = 'cm')


# Create the Spearman plot
spearman_plot <- ggplot(subset(long_results, variable == "Spearman"), 
                        aes(x = Tissue, y = value, fill = Tissue)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.4) +
  labs(title = "Spearman",
       y = "Spearman Correlation") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 1)) +
  scale_fill_manual(name = "Tissue", values = colors)

ggsave("./figures/Spearman_all.pdf", width = 10, height = 8, units = 'cm')
ggsave("./figures/Spearman_all.png", width = 10, height = 8, units = 'cm')


# Deviation of RMSE, Pearson and Spearman correlation for cell types -------------------
### USE NEXT PLOT TO SHOW THE DEVIATION ACROSS CELLTPES ###

# Read csv files containing RMSE, Pearson and Spearman correlations for overall tissue
results_kidney <- read.table("./local/kidney/evaluation_results_all_kidney.csv", header = TRUE)
results_lung <- read.table("./local/lung/evaluation_results_all_lung.csv", header = TRUE)
results_breast <- read.table("./local/breast/evaluation_results_all_breast.csv", header = TRUE)

# Read csv files containing RMSE, Pearson and Spearman correlations for each cell type
results_kidney_celltype <- read.csv("./local/kidney/evaluation_results_celltype_kidney.csv")
results_lung_celltype <- read.csv("./local/lung/evaluation_results_celltype_lung.csv")
results_breast_celltype <- read.csv("./local/breast/evaluation_results_celltype_breast.csv")

# Add a column to each dataset to identify the tissue
results_kidney$Tissue <- "Kidney"
results_lung$Tissue <- "Lung"
results_breast$Tissue <- "Breast"

results_kidney_celltype$Tissue <- "Kidney"
results_lung_celltype$Tissue <- "Lung"
results_breast_celltype$Tissue <- "Breast"

# Combine the overall and cell type datasets for all tissues
results_all <- rbind(results_kidney, results_lung, results_breast)
results_celltypes <- rbind(results_kidney_celltype, results_lung_celltype, results_breast_celltype)

# Calculate Absolute Deviations
results_celltypes_deviation <- results_celltypes %>%
  pivot_longer(cols = c(RMSE, Pearson, Spearman), names_to = "Metric", values_to = "Value") %>%
  left_join(results_all %>% select(Tissue, RMSE, Pearson, Spearman), by = "Tissue") %>%
  mutate(AbsoluteDeviation = case_when(
    Metric == "RMSE" ~ abs(Value - RMSE),
    Metric == "Pearson" ~ abs(Value - Pearson),
    Metric == "Spearman" ~ abs(Value - Spearman)
  ))

# Plot for Kidney
ggplot(data = filter(results_celltypes_deviation, Tissue == "Kidney"), aes(x = cell_type, y = AbsoluteDeviation, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Kidney", x = NULL, y = "Absolute Deviation from Average") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels
    legend.title = element_blank()) +
  scale_fill_brewer(palette = "Set1")

ggsave("./figures/deviation_celltypes_kidney.pdf", width = 12, height = 10, units = 'cm')
ggsave("./figures/deviation_celltypes_kidney.png", width = 12, height = 10, units = 'cm')


# Lung Plot
ggplot(data = filter(results_celltypes_deviation, Tissue == "Lung"), aes(x = cell_type, y = AbsoluteDeviation, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Lung", x = NULL, y = "Absolute Deviation from Average") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels
    legend.title = element_blank()  # Remove legend title
  ) +
  scale_fill_brewer(palette = "Set1")

ggsave("./figures/deviation_celltypes_lung.pdf", width = 12, height = 10, units = 'cm')
ggsave("./figures/deviation_celltypes_lung.png", width = 12, height = 10, units = 'cm')


# Breast Plot
ggplot(data = filter(results_celltypes_deviation, Tissue == "Breast"), aes(x = cell_type, y = AbsoluteDeviation, fill = Metric)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Breast", x = NULL, y = "Absolute Deviation from Average") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  # Center the title
    axis.text.x = element_text(angle = 90, hjust = 1),  # Rotate x-axis labels
    legend.title = element_blank()  # Remove legend title
  ) +
  scale_fill_brewer(palette = "Set1")

ggsave("./figures/deviation_celltypes_breast.pdf", width = 12, height = 10, units = 'cm')
ggsave("./figures/deviation_celltypes_breast.png", width = 12, height = 10, units = 'cm')


# RMSE, Pearson and Spearman for each cell type ---------------------------

# Read csv files containing RMSE, Pearson, and Spearman correlations
results_kidney <- read.table("./local/kidney/evaluation_results_all_kidney.csv", header = TRUE)
results_lung <- read.table("./local/lung/evaluation_results_all_lung.csv", header = TRUE)
results_breast <- read.table("./local/breast/evaluation_results_all_breast.csv", header = TRUE)

# Read cell type data
results_kidney_celltype <- read.csv("./local/kidney/evaluation_results_celltype_kidney.csv")
results_lung_celltype <- read.csv("./local/lung/evaluation_results_celltype_lung.csv")
results_breast_celltype <- read.csv("./local/breast/evaluation_results_celltype_breast.csv")

# Add a column to distinguish between overall and cell types
results_kidney$cell_type <- "Overall"
results_lung$cell_type <- "Overall"
results_breast$cell_type <- "Overall"

# Combine overall and cell type data for each tissue
results_kidney_combined <- rbind(results_kidney, results_kidney_celltype)
results_lung_combined <- rbind(results_lung, results_lung_celltype)
results_breast_combined <- rbind(results_breast, results_breast_celltype)

# Reshape data to long format for plotting
long_results_kidney <- melt(results_kidney_combined, id.vars = "cell_type", 
                            measure.vars = c("RMSE", "Pearson", "Spearman"))
long_results_lung <- melt(results_lung_combined, id.vars = "cell_type", 
                          measure.vars = c("RMSE", "Pearson", "Spearman"))
long_results_breast <- melt(results_breast_combined, id.vars = "cell_type", 
                            measure.vars = c("RMSE", "Pearson", "Spearman"))

# Define colors for each cell type in kidney
colors_kidney <- c(
  "B" = "#7B90AA",
  "Duct epithelial cell" = "#7D7D61",
  "Kidney epithelial cell" = "#D9C6C3",
  "Meso-epithelial cell" = "#B1BECD",
  "Monocyte" = "#C4B893",
  "Overall" = "#BF878A", 
  "Pericytes" = "#E6E4C5",
  "T" = "#BFBAE5"
)

# Plot cell types in kidney
ggplot(long_results_kidney, aes(x = variable, y = value, fill = cell_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.4) +
  labs(title = "Kidney", y = "", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(0, 1)) + 
  scale_fill_manual(name = "Cell Type", values = colors_kidney)

ggsave("./figures/kidney_celltypes.pdf", width = 15, height = 8, units = 'cm')
ggsave("./figures/kidney_celltypes.png", width = 15, height = 8, units = 'cm')

# Define colors for each cell type in lung
colors_lung <- c(
  "B" = "#7B90AA",
  "Fibroblast" = "#84A8A3",
  "Lung epithelial cell" = "#D9C6C3",
  "Mast cells" = "#82B9D7",
  "Meso-epithelial cell" = "#B1BECD",
  "Monocyte" = "#C4B893",
  "Overall" = "#BF878A",
  "T" = "#BFBAE5"
)

# Plot cell types in lung
ggplot(long_results_lung, aes(x = variable, y = value, fill = cell_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.4) +
  labs(title = "Lung", y = "", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(0, 1)) + 
  scale_fill_manual(name = "Cell Type", values = colors_lung)

ggsave("./figures/lung_celltypes.pdf", width = 15, height = 8, units = 'cm')
ggsave("./figures/lung_celltypes.png", width = 15, height = 8, units = 'cm')


# Define colors for each cell type in breast
colors_breast <- c(
  "Duct epithelial cell" = "#7D7D61",
  "Fibroblast" = "#84A8A3",
  "Meso-epithelial cell" = "#B1BECD",
  "Monocyte" = "#C4B893",
  "Overall" = "#BF878A"
)

# Plot cell types in breast
ggplot(long_results_breast, aes(x = variable, y = value, fill = cell_type)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.5), width = 0.4) +
  labs(title = "Breast", y = "", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.title.x = element_blank(), 
        legend.position = "right",
        plot.title = element_text(hjust = 0.5)) + 
  scale_y_continuous(limits = c(0, 1)) + 
  scale_fill_manual(name = "Cell Type", values = colors_breast)


ggsave("./figures/breast_celltypes.pdf", width = 15, height = 8, units = 'cm')
ggsave("./figures/breast_celltypes.png", width = 15, height = 8, units = 'cm')


# Stacked barplot of cell count in each tissue ----------------------------

# Prepare data with annotation level 1
kidney <- readRDS("./local/kidney/kidneyRep2.rds")
kidney <- prepare_data(kidney,
                       integrate_data = FALSE,
                       annotation_selfCluster = TRUE, 
                       remove_empty_droplets = FALSE,
                       anno_level = 1)
saveRDS(kidney ,"./local/kidney/kidneyRep2_prepared_anno1.rds")

breast <- readRDS("./local/breast/breast.rds")
breast <- prepare_data(breast,
                       integrate_data = FALSE,
                       annotation_selfCluster = TRUE, 
                       remove_empty_droplets = FALSE,
                       anno_level = 1)
saveRDS(breast ,"./local/breast/breast_prepared_anno1.rds")

lung <- readRDS("./local/lung/lung.rds")
lung <- prepare_data(lung,
                       integrate_data = FALSE,
                       annotation_selfCluster = TRUE, 
                       remove_empty_droplets = FALSE,
                       anno_level = 1)
saveRDS(lung ,"./local/lung/lung_prepared_anno1.rds")

# Read level 1 annotated files
kidney <- readRDS("./local/kidney/kidneyRep2_prepared_anno1.rds")
breast <- readRDS("./local/breast/breast_prepared_anno1.rds")
lung <- readRDS("./local/lung/lung_prepared_anno1.rds")

# Extract cell type information and tissue labels
get_celltype_data <- function(seurat_obj, tissue) {
  data.frame(
    cell_type = seurat_obj@meta.data$cell_type,
    tissue = tissue
  )
}

# Combine data from all tissues
kidney_data <- get_celltype_data(kidney, "Kidney")
breast_data <- get_celltype_data(breast, "Breast")
lung_data <- get_celltype_data(lung, "Lung")
celltype_data <- bind_rows(kidney_data, breast_data, lung_data)

# Calculate the percentage of each cell type within each tissue
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

ggsave("./figures/cell_count_stacked_barplot.pdf", width = 11, height = 10, units = 'cm')
ggsave("./figures/cell_count_stacked_barplot.png", width = 11, height = 10, units = 'cm')
