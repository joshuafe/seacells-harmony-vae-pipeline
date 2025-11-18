# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org/"))

# Install and load packages
if (!requireNamespace("flowCore", quietly = TRUE)) {
    BiocManager::install("flowCore")
}
if (!requireNamespace("uwot", quietly = TRUE)) {
    install.packages("uwot")
}
if (!requireNamespace("ggplot2", quietly = TRUE)) {
    install.packages("ggplot2")
}
if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
}
if (!requireNamespace("tidyr", quietly = TRUE)) {
    install.packages("tidyr")
}
if (!requireNamespace("leiden", quietly = TRUE)) {
    install.packages("leiden")
}
if (!requireNamespace("RANN", quietly = TRUE)) {
    install.packages("RANN")
}
if (!requireNamespace("igraph", quietly = TRUE)) {
    install.packages("igraph")
}

library(flowCore)
library(uwot)
library(ggplot2)
library(dplyr)
library(tidyr)
library(leiden)
library(RANN)
library(igraph)

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
fcs_file <- args[1]
output_dir <- args[2]

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load FCS file
fcs <- read.FCS(fcs_file, transformation = FALSE)

# Get expression data
exprs_data <- exprs(fcs)

# Define channels to exclude
channels_to_exclude <- c("FSC-A", "FSC-H", "SSC-A", "SSC-H", "Time", "Viability")
channels_to_transform <- colnames(exprs_data)[!colnames(exprs_data) %in% channels_to_exclude]

# Arcsinh transformation
asinh_trans <- arcsinhTransform(a = 0, b = 1/5)
for (channel in channels_to_transform) {
    if (channel %in% colnames(exprs_data)) {
        exprs_data[, channel] <- asinh_trans(exprs_data[, channel])
    }
}

# Remove zero-variance columns
variances <- apply(exprs_data, 2, var)
zero_var_cols <- names(variances[variances == 0])
if (length(zero_var_cols) > 0) {
    print(paste("Removing zero-variance columns:", paste(zero_var_cols, collapse = ", ")))
    exprs_data <- exprs_data[, !(colnames(exprs_data) %in% zero_var_cols)]
}

if (!requireNamespace("RANN", quietly = TRUE)) {
    install.packages("RANN")
}
library(RANN)

# ... (previous code)

# PCA
print("Running PCA...")
pca_out <- prcomp(exprs_data, center = TRUE, scale. = TRUE)
n_pcs <- min(30, ncol(pca_out$x))
pca_data <- pca_out$x[, 1:n_pcs]

if (!requireNamespace("igraph", quietly = TRUE)) {
    install.packages("igraph")
}
library(igraph)

# ... (previous code)

# Leiden clustering
print("Running Leiden clustering...")
nn_graph <- RANN::nn2(pca_data, k = 30)$nn.idx
# Create a graph from the nearest-neighbor list
adj <- matrix(0, nrow = nrow(pca_data), ncol = nrow(pca_data))
for (i in 1:nrow(pca_data)) {
    adj[i, nn_graph[i, -1]] <- 1
}
graph <- graph_from_adjacency_matrix(adj, mode = "undirected")
leiden_clusters <- leiden(graph, resolution_parameter = 0.1)
umap_df <- as.data.frame(pca_data)
umap_df$leiden_cluster <- as.factor(leiden_clusters)

# UMAP
print("Running UMAP on PCA data...")
set.seed(42)
umap_out <- umap(pca_data, n_neighbors = 15, min_dist = 0.1, n_components = 2)
umap_df$UMAP1 <- umap_out[,1]
umap_df$UMAP2 <- umap_out[,2]


# Plot UMAP colored by Leiden clusters
print("Generating UMAP plot colored by Leiden clusters...")
ggplot(umap_df, aes(x = UMAP1, y = UMAP2, color = leiden_cluster)) +
    geom_point(size = 1, alpha = 0.7) +
    labs(title = "UMAP colored by Leiden Clusters") +
    theme_minimal()
ggsave(file.path(output_dir, "umap_leiden_clusters.png"), width = 10, height = 8)


# Marker plots
markers <- c("CD3+CD19", "CD56", "CD16+TIGIT", "CD4+CD33", "CD8+CD14", "CD45")

for (marker in markers) {
    if (marker %in% colnames(exprs_data)) {
        print(paste("Generating plot for", marker))
        
        p <- ggplot(umap_df, aes_string(x = "UMAP1", y = "UMAP2", color = paste0("`", marker, "`"))) +
            geom_point(size = 1, alpha = 0.7) +
            scale_color_viridis_c() +
            labs(title = paste("UMAP colored by", marker, "expression")) +
            theme_minimal()
        
        ggsave(file.path(output_dir, paste0("umap_", gsub("[+]", "_", marker), ".png")), plot = p, width = 10, height = 8)
    }
}

print("Script finished.")