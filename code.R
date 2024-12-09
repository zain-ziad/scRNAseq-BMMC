# Assignment 01 : scRNAseq Worksheet
# Course : Single Cell Bioinformatics 2024-25
# Submitted by: Zain Ziad (7068669)

# setting seed for reproducibility
set.seed(11)

# loading necessary libraries
suppressPackageStartupMessages({
library(dplyr)
library(spatstat.core)
library(Seurat)
library(patchwork)
library(DoubletFinder)
library(SingleR)
library(enrichR)
library(SingleCellExperiment)
library(SeuratWrappers)
library(tidyverse)
library(celldex)
library(EnhancedVolcano)
library(monocle3)
library(CellChat)
})

# make figure directories
dir.create("figures")
dir.create("figures/After")

### TASK 1: loading the data

# define directory containing BMMC and CD34 data
directory <- "scbi_ds1//"

# define the file names to be loaded
file_names <- c("GSM4138872_scRNA_BMMC_D1T1.rds",
                "GSM4138873_scRNA_BMMC_D1T2.rds",
                "GSM4138874_scRNA_CD34_D2T1.rds",
                "GSM4138875_scRNA_CD34_D3T1.rds")


load_dataset <- function(directory, file_name) {
  # a func for loading a sample and turning it into a seurat Object

  # construct the directory of the file to be loaded
  file_directory <- paste(directory, file_name, sep = "")

  # read the RDS file containing the raw counts
  raw_counts <- readRDS(file = file_directory)

  # create a seurat object from the raw counts
  seurat_obj <- CreateSeuratObject(
                counts = raw_counts,
                project = file_name,
                assay = "RNA"
            )
  return(seurat_obj)
}

# init an empty list to store seurat objs
list_of_samples <- list()

# use load_dataset on each of the samples
for (i in 1:length(file_names)) {
  seurat_obj <- load_dataset(directory, file_names[i])
  list_of_samples[[i]] <- seurat_obj
}

### TASK 2: create the sample sheet

sample_sheet <- data.frame(
  file = file_names,
  names = c("BMMC_D1T1", "BMMC_D1T2", "CD34_D2T1", "CD34_D3T1"),
  donor = c("D1", "D1", "D2", "D3"),
  replicate = c("T1", "T2", "T1", "T1"),
  sex = c("F", "F", "M", "F"),
  group = c("BMMC", "BMMC", "CD34", "CD34")
)

add_metadata <- function(seurat_obj, sample_sheet, i) {
  # a function to add all of the metadata from
  # the sample_sheet to the Seurat obj

  seurat_obj$orig.ident <- sample_sheet$names[i]
  seurat_obj$donor <- sample_sheet$donor[i]
  seurat_obj$replicate <- sample_sheet$replicate[i]
  seurat_obj$sex <- sample_sheet$sex[i]
  seurat_obj$group <- sample_sheet$group[i]
  return(seurat_obj)
}

### TASK 3: add Meta-data

# apply metadata to all samples
list_of_samples <- lapply(1:length(list_of_samples), function(i) {
  add_metadata(list_of_samples[[i]], sample_sheet, i)
})

# to answer the asked questions
for (sample in list_of_samples) {
    print(seurat_obj@project.name)

    print("How many cells are in this sample?")
    print(ncol(sample))

    print("How many genes are in this expression matrice?")
    print(nrow(sample))

    print("What information is now part of the meta-data of the object?")
    print(head(sample@meta.data))
}

### TASK 4: preprocessing

make_qcplots <- function(sample, save_dir, calc_pctmt = TRUE) {

    # a function that makes the basic qc
    # violin and scatter plots for each of the sample
    # calc_pctmt controls either to calulate mitochondrial
    # percentage or not

    sample_name <- sample@project.name

    if (calc_pctmt) {
        features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
    } else {
        features <- c("nFeature_RNA", "nCount_RNA")
    }

    qc_plt <- VlnPlot(sample,
    features = features,
    ncol = 3)

    qc_scatter <- FeatureScatter(sample, feature1 = "nCount_RNA",
                                feature2 = "nFeature_RNA")

    output_plt <- qc_plt / qc_scatter

    sample_name <- sample@project.name
    ggsave(filename = paste0(save_dir, "/QC_", sample_name, ".png"),
    plot = output_plt, width = 15, height = 12, dpi = 100, bg = "white")
}

for (sample in list_of_samples) {
    # calculate MT percentage
    sample[["percent.mt"]] <- PercentageFeatureSet(sample, pattern = "^MT-")

    make_qcplots(sample, "figures")
    }

# filter by qc metrics
for (i in 1:length(list_of_samples)) {

    list_of_samples[[i]] <- subset(list_of_samples[[i]],
        subset = nFeature_RNA > 750 &
                    nFeature_RNA < 2500 &
                    nCount_RNA > 1500 &
                    nCount_RNA < 6000)

    make_qcplots(list_of_samples[[i]], "figures/After", FALSE)
    }

# normalize each sample
for (i in 1:length(list_of_samples)) {
    list_of_samples[[i]] <- NormalizeData(list_of_samples[[i]])
}

# doublet prediction using DoubletFinder
doublet_pred <- function(sample) {

    # a function that runs the doubletfinder pipline
    # on each of the samples
    # returns the processed seurat obj along with
    # making barplots for the distribution of doublets
    # vs singlets in the dataset

    sample_name <- sample@project.name

    # feature selection
    sample <- FindVariableFeatures(sample)
    sample <- ScaleData(sample)
    sample <- RunPCA(sample)

    # run parameter sweep
    sweep.res <- paramSweep(sample, PCs = 1:10)

    # summarize results
    sweep.stats <- summarizeSweep(sweep.res)

    # find best pK value
    bcmvn <- find.pK(sweep.stats)
    max_bcmetric_index <- which.max(bcmvn$BCmetric)
    pk_val <- bcmvn$pK[max_bcmetric_index]
    print("Using pK value: ")
    print(pk_val)

    # calculate expected doublets using poisson estimation
    nexp_poi <- round(0.075 * nrow(sample@meta.data))

    # run doubletfinder
    sample <- doubletFinder(sample,
                            PCs = 1:10, pN = 0.25,
                            pK = as.numeric(as.character(pk_val)),
                            nExp = nexp_poi, reuse.pANN = FALSE, sct = FALSE)

    # next we will visualize it
    # get the DoubletFinder results column name
    df_name <- colnames(sample@meta.data)[grep("DF.classification",
                                            colnames(sample@meta.data))]

    # calculate counts
    doublet_counts <- table(sample@meta.data[,df_name])

    # convert to data frame for ggplot
    df_counts <- data.frame(
    Status = names(doublet_counts),
    Count = as.numeric(doublet_counts)
    )

    # barplot
    p1 <- ggplot(df_counts, aes(x = Status, y = Count, fill = Status)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = Count), vjust = -0.5) +
    theme_minimal() +
    labs(title = "Distribution of Singlets vs Doublets",
        y = "Number of Cells",
        x = "") +
    scale_fill_manual(values = c("Doublet" = "#FC8D62",
                                 "Singlet" = "#66C2A5")) +
    theme(legend.position = "none")

    ggsave(filename = paste0("figures", "/DFResults_", sample_name, ".png"),
    plot = p1, width = 6, height = 5, dpi = 100, bg = "white")

    sample <- subset(sample,
                cells = rownames(sample@meta.data)[sample@meta.data[,df_name] == "Singlet"]) # nolint

    return(sample)
}

# run the function on each samples
for (i in 1:length(list_of_samples)) {
    list_of_samples[[i]] <- doublet_pred(list_of_samples[[i]])
}

### TASK 4.2, 5.1, 5.2:
### batch correction, dimensionality reduction, clustering

# (1)

pca_clusters <- function(sample, find_clusters) {
    # a function that runs the pca and clustring (if find_clusters),
    # returns the processed seurat object
    # in addtion to creating elbow plots for pca

    sample <- FindVariableFeatures(sample)
    sample <- ScaleData(sample)
    sample <- RunPCA(sample)

    sample_elbow <- ElbowPlot(sample, ndims = 50, reduction = "pca")
    ggsave(filename = paste0("figures/", sample@project.name, "_elbow", ".png"),
    plot = sample_elbow, width = 5, height = 5, dpi = 100, bg = "white")

    sample <- FindNeighbors(sample, dims = 1:30, reduction = "pca")

    if (find_clusters) {
        sample <- FindClusters(sample, resolution = 0.2)
    }

    # umap for visualization
    sample <- RunUMAP(sample, dims = 1:30, reduction = "pca")

    if (find_clusters) {
        cluster_map <- DimPlot(sample,
                                reduction = "umap",
                                group.by = "seurat_clusters",
                                label = TRUE)
        ggsave(filename = paste0("figures/", sample@project.name,
                        "_clustermap", ".png"),
                        plot = cluster_map, width = 7, height = 5, dpi = 100, bg = "white")
    }

    return(sample)
}

# merge all samples without integration

merged_seurat <- merge(list_of_samples[[1]],
                      y = list_of_samples[2:4],
                      add.cell.ids = names(list_of_samples),
                      project = "merged_project")

# run the clustering steps on merged dataset
merged_seurat <- pca_clusters(merged_seurat, find_clusters = FALSE)

# (2)

# seurat v4 workflow... find integration anchors
anchors <- FindIntegrationAnchors(list_of_samples,
                                 scale = FALSE)

# integrate data
integrated_seurat <- IntegrateData(anchorset = anchors)
integrated_seurat@project.name <- "integrated_project"

# now we switch to integrated assay for downstream analysis
DefaultAssay(integrated_seurat) <- "integrated"

# do QC again
make_qcplots(integrated_seurat, "figures")
# no further qc was required

# run the clustering steps on integrated dataset
integrated_seurat <- pca_clusters(integrated_seurat, find_clusters = TRUE)

# to compare batch effects
p1 <- DimPlot(merged_seurat,
              reduction = "umap",
              group.by = "donor") + ggtitle("Merged")

p2 <- DimPlot(integrated_seurat,
              reduction = "umap",
              group.by = "donor") + ggtitle("Integrated")

ggsave(filename = paste0("figures/", "compare_batcheffects.png"),
        plot = (p1 + p2), width = 12, height = 5, dpi = 100, bg = "white")

### TASK 6:

# 6.1: automatic cell type annotation
reference_data <- HumanPrimaryCellAtlasData()

# since singler expects the expression data
# get expression matrix from seurat object
singler_data <- GetAssayData(integrated_seurat, slot = "data")

predictions <- SingleR(test = singler_data, ref = reference_data,
                       labels = reference_data$label.main)

# now we can add these predictions to seurat object
integrated_seurat$singler.labels <- predictions$labels

# extract some plots
singler_clusters <- DimPlot(integrated_seurat,
                    reduction = "umap",
                    group.by = "singler.labels")

ggsave(filename = paste0("figures", "/singler_clusters.png"),
plot = singler_clusters, width = 7, height = 5, dpi = 100, bg = "white")

singler_score <- plotScoreHeatmap(predictions)
ggsave(filename = paste0("figures", "/singler_scores.png"),
plot = singler_score, width = 10, height = 7, dpi = 100, bg = "white")

# 6.2: manual annotation
celltype_markers <- FindAllMarkers(object = integrated_seurat)

# groupby clusters and keep only the top 25 genes
# for each cluster with significant pval

celltype_markers <- group_by(celltype_markers, cluster)
celltype_markers <- filter(celltype_markers, p_val_adj < 0.05)
celltype_markers <- slice_max(celltype_markers, n = 25, order_by = avg_log2FC)

# visualization with a heatmap
top3 <- slice_max(celltype_markers, n = 3, order_by = avg_log2FC)
top3 <- ungroup(top3)

gene_rankings <- DoHeatmap(integrated_seurat, features = top3$gene) + NoLegend()
ggsave(filename = paste0("figures", "/gene_rankings.png"),
plot = gene_rankings, width = 15, height = 15, dpi = 100, bg = "white")

# now we look at these genes closely
# umap and violin for top 5 cluster genes

for (c in unique(top3$cluster)) {
    # get genes for this specific cluster
    cluster_genes <- filter(top3, cluster == c)
    cluster_genes <- pull(cluster_genes, gene)

    p1 <- FeaturePlot(integrated_seurat,
                    features = cluster_genes)
    p2 <- VlnPlot(integrated_seurat,
                    features = cluster_genes)

    ggsave(filename = paste0("figures", "/clusterplot_", c, ".png"),
    plot = p1, width = 12, height = 8, dpi = 100, bg = "white")

    ggsave(filename = paste0("figures", "/vlnplot_", c, ".png"),
    plot = p2, width = 12, height = 6, dpi = 100, bg = "white")
}

# using the reference cell marker table
# we can annotate the cell types to the clusters like so

integrated_seurat$cell_types <- plyr::mapvalues(
  integrated_seurat$seurat_clusters,
  from = c(0:13),
  to = c("CD14+ Monocytes",
         "T Cells",
         "GMP/Neutrophils",
         "Common Lymphoid Progenitor",
         "Common Lymphoid Progenitor",
         "Erythrocytes",
         "Pre B cells Progenitor",
         "CD8+ T Cells",
         "NK cells",
         "plasmacytoid Dendritic Cells",
         "B Cells",
         "Pre B cells Progenitor",
         "NK cells",
         "Monocytes")
)

celltype_clusters <- DimPlot(integrated_seurat,
                    reduction = "umap",
                    group.by = "cell_types")

ggsave(filename = paste0("figures", "/celltype_clusters.png"),
plot = celltype_clusters, width = 13, height = 8, dpi = 100, bg = "white")

# 6.3: cell type proportions
integrated_seurat$sample <- paste0(integrated_seurat$donor,
                                    integrated_seurat$replicate)

# calculate proportions for each cell type per sample
proportions <- table(integrated_seurat$sample, integrated_seurat$cell_types)
proportions <- proportions / rowSums(proportions)
proportions <- as.data.frame(proportions)
colnames(proportions) <- c("Sample", "CellType", "Proportion")


# making a stacked bar plot
p1 <- ggplot(proportions, aes(x = Sample, y = Proportion, fill = CellType)) +
    geom_bar(stat = "identity", position = "fill") +
    theme_minimal() +
    labs(x = "Sample",
        y = "Proportion",
        title = "Cell Type Proportions by Sample")

ggsave(filename = paste0("figures", "/celltype_proportions.png"),
plot = p1, width = 10, height = 8, dpi = 100, bg = "white")

# (7.1) differential expression

bcell_tcell <- FindMarkers(integrated_seurat,
                             ident.1 = "B Cells",
                             ident.2 = c("T Cells", "CD8+ T Cells"),
                             group.by = "cell_types")

bcell_tcell <- subset(bcell_tcell, p_val_adj < 0.05)

bcell_tcell_p <-  EnhancedVolcano(bcell_tcell,
    lab = rownames(bcell_tcell),
    x = 'avg_log2FC',
    y = 'p_val_adj',
    title = 'B Cells vs T Cells',
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize = 6.0)

ggsave(filename = paste0("figures", "/bcellvtcell.png"),
plot = bcell_tcell_p, width = 10, height = 8, dpi = 100, bg = "white")


# tcells vs monocytes

tcell_mono <- FindMarkers(integrated_seurat,
                             ident.1 = c("Monocytes", "CD14+ Monocytes"),
                             ident.2 = c("T Cells", "CD8+ T Cells"),
                             group.by = "cell_types")

tcell_mono <- subset(tcell_mono, p_val_adj < 0.05)

tcell_mono_p <-  EnhancedVolcano(tcell_mono,
    lab = rownames(tcell_mono),
    x = 'avg_log2FC',
    y = 'p_val_adj',
    title = 'Monocytes vs T Cells',
    FCcutoff = 0.5,
    pointSize = 3.0,
    labSize = 6.0)

ggsave(filename = paste0("figures", "/tcellvmono.png"),
plot = tcell_mono_p, width = 10, height = 8, dpi = 100, bg = "white")

# (7.2) plot differentially expressed genes

Idents(integrated_seurat) <- "cell_types"
de_genes <- FindAllMarkers(integrated_seurat)

# transform p_val_adj and handle zeros
small_value <- 1e-300
de_genes <- de_genes %>%
  mutate(
    p_val_adj = ifelse(p_val_adj == 0, small_value, p_val_adj),
    significance = -log10(p_val_adj) #so less = more
  )

# we get the top5 genes in each celltype acc to log2FC
de_genes_top5 <- de_genes %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5) %>%
  ungroup()

# then we filter to just those genes
filtered_de_genes <- filter(de_genes, gene %in% de_genes_top5$gene)

# dot plot
p1 <- ggplot(filtered_de_genes, aes(x = cluster, y = gene)) +
  geom_point(aes(size = significance, color = avg_log2FC)) +
  scale_color_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "Avg Log2FC"
  ) +
  scale_size_continuous(
    range = c(3, 10),  # Adjust sizes as needed
    name = "-log10(p_adj)"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Dotplot of Top 5 Genes per Cluster",
    x = "Cell Type",
    y = "Genes"
  )

ggsave(filename = paste0("figures", "/degenes_per_celltype.png"),
plot = p1, width = 10, height = 12, dpi = 100, bg = "white")

# (8) pathway analysis

# 8.1
Idents(integrated_seurat) <- "group"

bmmc_cd34 <- FindMarkers(integrated_seurat,
                             ident.1 = "BMMC",
                             ident.2 = "CD34",
                             group.by = "group")

bmmc_cd34 <- subset(bmmc_cd34, p_val_adj < 0.05)
bmmc_cd34_f <- slice_max(bmmc_cd34, n = 5, order_by = avg_log2FC)

write.csv(bmmc_cd34_f, "bmmc_cd34.csv", row.names = TRUE)

# for only monocyte cells

monocytes <- subset(integrated_seurat,
                    subset = cell_types %in% c("CD14+ Monocytes", "Monocytes"))

bmmc_cd34_monocytes <- FindMarkers(
  monocytes,
  ident.1 = "BMMC",
  ident.2 = "CD34",
  group.by = "group"
)

bmmc_cd34_monocytes <- subset(bmmc_cd34_monocytes, p_val_adj < 0.05)
bmmc_cd34_monocytes <- slice_max(bmmc_cd34_monocytes,
                        n = 5, order_by = avg_log2FC)

write.csv(bmmc_cd34_monocytes, "bmmc_cd34_monocytes.csv", row.names = TRUE)

# (8.2) pathway analysis

GO_Biological_Process <- DEenrichRPlot(
  integrated_seurat,
  ident.1 = "BMMC",
  ident.2 = "CD34",
  p.val.cutoff = 0.05,
  enrich.database = "GO_Biological_Process_2021",
  numPathways = 10,
  max.genes = 200
)

ggsave(filename = paste0("figures", "/GO_Biological_Process.png"),
plot = GO_Biological_Process, width = 12, height = 15, dpi = 100, bg = "white")

### TASK 9: trajectory analysis

#since we are interested in bcell trajectory
btraj_seurat <- subset(integrated_seurat,
                    cell_types %in% c("Common Lymphoid Progenitor",
                                    "Pre B cells Progenitor",
                                    "B Cells"))

# setting up the data in monocle3
bcell_cds <- as.cell_data_set(btraj_seurat)
bcell_cds <- cluster_cells(bcell_cds)
bcell_cds <- learn_graph(bcell_cds, use_partition = FALSE)

# selecting root nodes manually
bcell_cds <- order_cells(bcell_cds)

# plot trajectories based on celltype
p1 <- plot_cells(
  cds = bcell_cds,
  color_cells_by = "cell_types",
  show_trajectory_graph = FALSE,
  group_label_size = 4
)

ggsave(filename = paste0("figures/", "bcell_trajcelltype.png"),
    plot = p1, width = 8, height = 10, dpi = 300, bg = "white")

# color based on pseudotime
p1 <- plot_cells(
  cds = bcell_cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

ggsave(filename = paste0("figures/", "bcell_trajptime.png"),
    plot = p1, width = 8, height = 10, dpi = 300, bg = "white")

# selecting root nodes automatically

# using monocle3 helper function to find nodes
# automatically

get_earliest_principal_node <- function(cds, start_cell_type) {
  cell_ids <- which(colData(cds)$cell_types == start_cell_type)

  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]

  root_pr_nodes
}

bcell_cds <- order_cells(bcell_cds,
            root_pr_nodes = get_earliest_principal_node(bcell_cds, "Common Lymphoid Progenitor"))

p1 <- plot_cells(
  cds = bcell_cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE
)

ggsave(filename = paste0("figures/", "auto_bcell_trajptime.png"),
    plot = p1, width = 8, height = 10, dpi = 300, bg = "white")

### TASK 10:  cell-cell communication
common_celltypes <- intersect(unique(integrated_seurat$cell_types[integrated_seurat$group == "CD34"]),
                            unique(integrated_seurat$cell_types[integrated_seurat$group == "BMMC"]))

common_seurat <- subset(integrated_seurat, cell_types %in% common_celltypes)

cellchat_process <- function(cc_object) {
  # select cellchat database
  cc_db <- CellChatDB.human

  # preprocess
  cc_object@DB <- cc_db
  cc_object <- subsetData(cc_object)
  future::plan("multisession", workers = 12)
  cc_object <- identifyOverExpressedGenes(cc_object)
  cc_object <- identifyOverExpressedInteractions(cc_object)
  cc_object@idents <- droplevels(cc_object@idents)

  # process
  cc_object <- computeCommunProb(cc_object)
  cc_object <- filterCommunication(cc_object, min.cells = 10)
  cc_object <- computeCommunProbPathway(cc_object)
  cc_object <- aggregateNet(cc_object)

  return(cc_object)
}

# create seurat obj for each group
bmmc_data <- subset(common_seurat, group == "BMMC")
cd34_data <- subset(common_seurat, group == "CD34")

# create cellchat obj
cc_bmmc <- createCellChat(object = bmmc_data,
                            group.by = "cell_types", assay = "RNA")
                            
cc_cd34 <- createCellChat(object = cd34_data,
                            group.by = "cell_types", assay = "RNA")

# run the function on each group
cc_bmmc <- cellchat_process(cc_bmmc)
cc_cd34 <- cellchat_process(cc_cd34)

# find common pathways
pathways_bmmc <- cc_bmmc@netP$pathways
pathways_cd34 <- cc_cd34@netP$pathways

common_pathways <- intersect(pathways_bmmc, pathways_cd34)

# common pathways in both groups
print(common_pathways)

interaction_no_str <- function(cc_obj, name) {
  groupSize <- as.numeric(table(cc_obj@idents))
  png(paste0("figures/", name, "_interaction_number.png"),
        width = 800, height = 800, res = 100)
  num_inter <- netVisual_circle(cc_obj@net$count, vertex.weight = groupSize,
                  weight.scale = T, label.edge = F,
                  title.name = "Number of interactions") 
  dev.off()
  png(paste0("figures/", name, "_interaction_strength.png"),
      width = 800, height = 800, res = 100)
  str_inter <- netVisual_circle(cc_obj@net$weight, vertex.weight = groupSize,
                  weight.scale = T, label.edge = F,
                  title.name = "Interaction weights/strength")
  dev.off()
}

# make plots for
# number of interactions and
# interaction strength / weights
# for each group

interaction_no_str(cc_bmmc, "BMMC")
interaction_no_str(cc_cd34, "CD34")

# taking a closer look at MIF pathway
pathways.show <- "MIF"

png(paste0("figures/", "MIF_circleplot.png"),
        width = 800, height = 800, res = 100)

netVisual_aggregate(cc_bmmc,
                   signaling = pathways.show,
                   layout = "circle",
                   edge.weight.max = NULL,
                   edge.width.max = 8,
                   signaling.name = pathways.show)
dev.off()

# to make circle plots per group
groupSize <- as.numeric(table(cc_bmmc@idents))
mat <- cc_bmmc@netP$prob[,,"MIF"]

png(paste0("figures/", "MIF_per_group.png"),
        width = 2048, height = 2048, res = 100)
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

dev.off()