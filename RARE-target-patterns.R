#install.packages("Seurat")
# library(ggplot2)
# library(tidyr)
# library(dplyr)
library(Seurat)
library(tidyverse)

## LOADING DATA ----------------------
setwd("C:\\Users\\ayahs\\Downloads")
path.to.seurat.object <- "eye.rds"
eye.cluster <- read.csv("eye_clusters.csv")

eye.cluster <- eye.cluster %>% select(c("clust", "identity.super"))

# Load the Daniocell Seurat object and the latest annotations
daniocell <- readRDS(path.to.seurat.object)

## SUBSETTING DATA --------------------
# A function to re-generate the 'tissue subset' objects
# Usage: daniocell.tissue.subset(daniocell, tissue)
#    e.g. to recreate the axial mesoderm object: daniocell.tissue.subset(daniocell, "axial")
# Valid tissues are: blastomeres, periderm, axial, neural, hematopoietic, muscle, endoderm, pgc, eye, pronephros, 
#    glial, mesenchyme, epidermis, ionocytes, taste, otic, pigment, fin, mural
daniocell.tissue.subset <- function(daniocell, tissue) {
  if (!(tissue %in% names(daniocell@misc$tissue.umaps))) {
    message("The valid parameters for tissue are: ", paste0(names(daniocell@misc$tissue.umaps), collapse = ", "))
    stop("Please supply a valid tissue name.")
  }
  daniocell@reductions$umap@cell.embeddings <- as.matrix(daniocell@misc$tissue.umaps[[tissue]])
  cells.to.keep <- rownames(daniocell@meta.data)[daniocell@meta.data$subset.full == tissue]
  return(
    subset(daniocell, cells = cells.to.keep)
  )
}

# Get eye subset out from daniocell
daniocell.eye <- daniocell.tissue.subset(daniocell, "eye")

## PLOTTING CLUSTERINGS ---------------

# See what clusterings are available to plot
colnames(daniocell.eye@meta.data)

# Plot clusters
# DimPlot(daniocell.eye, group.by = "cluster", label = T)
# DimPlot(daniocell.eye, group.by = "hpf.nice", label = T)
# DimPlot(daniocell.eye, group.by = "stage.group", label = T)
# DimPlot(daniocell.eye, group.by = "identity.super.short", label = T)
# DimPlot(daniocell.eye, group.by = "identity.sub.short") # label = T doesn't work -- is it because there are empty values?

# Changing the active clustering
# Some plots don't let you configure which clustering to use, and will use the active one, which is changed with
#    the Idents command.
# Idents(daniocell.eye) <- "identity.super.short"

## TIME & ON/OFF --------------------------------
colnames(daniocell.eye@meta.data)
DimPlot(daniocell.eye, group.by = "hpf.nice", label = T)
# >=.3 --> "real"

# Generate list of clusters
cell.clusters <- FetchData(daniocell.eye, vars = "cluster")
cell.clusters <- tibble::rownames_to_column(cell.clusters)
colnames(cell.clusters) <- c("cell", "clust")
eye.clust.names <- left_join(cell.clusters, eye.cluster, by = "clust")

# on.col15a1b <- WhichCells(daniocell.eye, expression = col15a1b > 0.3)
# on.rx2 <- WhichCells(daniocell.eye, expression = rx2 > 0.3)
# off.ccnd1 <- WhichCells(daniocell.eye, expression = ccnd1 < 0.3)
after.12 <- WhichCells(daniocell.eye, expression = hpf >= 12)
after.14 <- WhichCells(daniocell.eye, expression = hpf >= 14)
after.24 <- WhichCells(daniocell.eye, expression = hpf >= 24)
after.36 <- WhichCells(daniocell.eye, expression = hpf >= 36)
after.48 <- WhichCells(daniocell.eye, expression = hpf >= 48)
after.60 <- WhichCells(daniocell.eye, expression = hpf >= 60)
after.72 <- WhichCells(daniocell.eye, expression = hpf >= 72)
after.84 <- WhichCells(daniocell.eye, expression = hpf >= 84)
after.96 <- WhichCells(daniocell.eye, expression = hpf >= 96)
bin.9 <- WhichCells(daniocell.eye, expression = hpf >= 108)
bin.8 <- setdiff(after.96, bin.9)
bin.7 <- setdiff(after.84, after.96)
bin.6 <- setdiff(after.72, after.84)
bin.5 <- setdiff(after.60, after.72)
bin.4 <- setdiff(after.48, after.60)
bin.3 <- setdiff(after.36, after.48)
bin.2 <- setdiff(after.24, after.36)
bin.1 <- setdiff(after.12, after.24)
#cells.cmz <- intersect(intersect(intersect(on.col15a1b, on.rx2), off.ccnd1), after.48)
#daniocell.eye@meta.data[cells.cmz, "ciliary.DE"] <- "cmz"
# creating the bins for hpf from the daniocell data

## GET ALL EYE GENES
daniocell.eye.genes <- rownames(GetAssayData(daniocell.eye, layer = "data"))
head(daniocell.eye.genes)

after.48 <- WhichCells(daniocell.eye, expression = hpf.nice >= 12)

## FUNCTION FOR SUBSETTING DANIOCELL BY CELL TYPE/CLUSTER
subsetCells <- function(daniocell.eye, clust.name) {
  # Reading in clusters file
  eye.cluster <- read.csv("eye_clusters.csv")
  eye.cluster <- eye.cluster %>% select(c("clust", "identity.super"))
  # Getting cluster names from daniocell
  cell.clusters <- FetchData(daniocell.eye, vars = "cluster")
  cell.clusters <- tibble::rownames_to_column(cell.clusters)
  colnames(cell.clusters) <- c("cell", "clust")
  eye.clust.names <- left_join(cell.clusters, eye.cluster, by = "clust")
  # Filter by given cluster name
  filtered.clust <- eye.clust.names %>% filter(clust == clust.name)
  # Get unique cell values in given cluster
  # filtered.cells <- unique(filtered.clust$cell)
  # print(filter.clust$cell)
  # Filter Daniocell object by given cluster
  daniocell.filtered <- subset(daniocell.eye, cells = filtered.clust$cell)
  
  return(
    daniocell.filtered
  )
}


## FUNCTION FOR DATASET GENERATION
# Make sure to run lines 70-87 before generateData()
generateData <- function(geneList, daniocell.filtered) {
  final_dataset <- tibble()

  for (selectedGene in geneList) {
    data_binned <- tibble()
    # append Bin1 data
      
      # append Bin10 data
      Bin1 <- FetchData(daniocell.filtered, vars = selectedGene, cells=bin.1)
      Bin1 <- tibble::rownames_to_column(Bin1)
      colnames(Bin1) <- c("cell", "geneExpression")
      Bin1 <- Bin1 %>% mutate(gene = selectedGene, bin = "bin.1")
    
      # append Bin2 data
      Bin2<- FetchData(daniocell.filtered, vars = selectedGene, cells=bin.2)
      Bin2 <- tibble::rownames_to_column(Bin2)
      colnames(Bin2) <- c("cell", "geneExpression")
      Bin2 <- Bin2 %>% mutate(gene = selectedGene, bin = "bin.2")
      
      # append Bin3 data
      Bin3<- FetchData(daniocell.filtered, vars = selectedGene, cells=bin.3)
      Bin3 <- tibble::rownames_to_column(Bin3)
      colnames(Bin3) <- c("cell", "geneExpression")
      Bin3 <- Bin3 %>% mutate(gene = selectedGene, bin = "bin.3")
      
      # append Bin4 data
      Bin4<- FetchData(daniocell.filtered, vars = selectedGene, cells=bin.4)
      Bin4 <- tibble::rownames_to_column(Bin4)
      colnames(Bin4) <- c("cell", "geneExpression")
      Bin4 <- Bin4 %>% mutate(gene = selectedGene, bin = "bin.4")
      
      # append Bin5 data
      Bin5 <- FetchData(daniocell.filtered, vars = selectedGene, cells=bin.5)
      Bin5 <- tibble::rownames_to_column(Bin5)
      colnames(Bin5) <- c("cell", "geneExpression")
      Bin5 <- Bin5 %>% mutate(gene = selectedGene, bin = "bin.5")
      
      # append Bin6 data
      Bin6 <- FetchData(daniocell.filtered, vars = selectedGene, cells=bin.6)
      Bin6 <- tibble::rownames_to_column(Bin6)
      colnames(Bin6) <- c("cell", "geneExpression")
      Bin6 <- Bin6 %>% mutate(gene = selectedGene, bin = "bin.6")
      
      # append Bin7 data
      Bin7 <- FetchData(daniocell.filtered, vars = selectedGene, cells=bin.7)
      Bin7 <- tibble::rownames_to_column(Bin7)
      colnames(Bin7) <- c("cell", "geneExpression")
      Bin7 <- Bin7 %>% mutate(gene = selectedGene, bin = "bin.7")
      
      # append Bin8 data
      Bin8 <- FetchData(daniocell.filtered, vars = selectedGene, cells=bin.8)
      Bin8 <- tibble::rownames_to_column(Bin8)
      colnames(Bin8) <- c("cell", "geneExpression")
      Bin8 <- Bin8 %>% mutate(gene = selectedGene, bin = "bin.8")
      
      # append Bin9 data
      Bin9 <- FetchData(daniocell.filtered, vars = selectedGene, cells=bin.9)
      Bin9 <- tibble::rownames_to_column(Bin9)
      colnames(Bin9) <- c("cell", "geneExpression")
      Bin9 <- Bin9 %>% mutate(gene = selectedGene, bin = "bin.9")
      
      # Combine all bins into one dataset with bin labels
      data_binned <- rbind(Bin1, Bin2) %>%
        rbind((.), Bin3) %>%
        rbind((.), Bin4) %>%
        rbind((.), Bin5) %>%
        rbind((.), Bin6) %>%
        rbind((.), Bin7) %>%
        rbind((.), Bin8) %>%
        rbind((.), Bin9)
  
      
      # Filter by "on" genes (expression > 0.3)
      #data_binned <- data_binned %>% filter(geneExpression > 0.3)
      # Refactor for data viz: put bins in proper order
      data_binned$bin <- factor(data_binned$bin, levels=c("bin.1", "bin.2", "bin.3", "bin.4", "bin.5", "bin.6", "bin.7", "bin.8", "bin.9"))
      
      # append data from new gene to previous dataset
      final_dataset <- rbind(final_dataset, data_binned)
    
  }
  return(
    final_dataset
  )
}

## GET GENES OF INTEREST-----------------------
genematrix <- scan("known_gene_names.txt", "")
intersected.genes <- intersect(genematrix, daniocell.eye.genes)

# test1 <- generateData(intersected.genes[1:20])
# first20 <- inner_join(test1, eye.clust.names, by = "cell")
# optic.cup <- first20 %>% filter(identity.super == "optic cup")
# optic.cells <- unique(optic.cup$cell)

daniocell.optic <- subsetCells(daniocell.eye, "eye.37")
daniocell.progens.18 <- subsetCells(daniocell.eye, "eye.18")

# daniocell.progens.2 <- subsetCells(daniocell.eye, "eye.2")
# daniocell.progens.26 <- subsetCells(daniocell.eye, "eye.26")
# daniocell.progens.11 <- subsetCells(daniocell.eye, "eye.11")

# LOOPING GENES ON INTEREST------------------
# Get list of "on" genes
filterAvgGeneExpression <- function(chunked.genes, daniocell.celltype) {
  # only include genes with mean gene expression > 0.3 in at least one bin
  on.genes <- c()
  for (var in chunked.genes) {
    print(var)
    filteredData <- generateData(var, daniocell.celltype)
    
    mean1 <- filteredData %>% 
      group_by(gene, bin) %>% 
      summarise(mean = mean(geneExpression))
    
    mean1 <- mean1 %>% filter(mean > 0.3)
    print(unique(mean1$gene))
    on.genes <- c(on.genes, unique(mean1$gene))
    print(on.genes)
  }
  return(on.genes)
}

filterMedianGeneExpression <- function(chunked.genes, daniocell.celltype) {
  # # Chunk the genes for functionality
  # chunked.genes <- split(geneList, ceiling(seq_along(intersected.genes)/20))
  # only include genes with median gene expression > 0.3 in at least one bin
  on.genes <- c()
  for (var in chunked.genes) {
    print(var)
    filteredData <- generateData(var, daniocell.celltype)
    print(head(filteredData))
    
    median1 <- filteredData %>% 
      group_by(gene, bin) %>% 
      summarise(med = median(geneExpression))
    
    median1 <- median1 %>% filter(med > 0.3)
    print(unique(median1$gene))
    on.genes <- c(on.genes, unique(median1$gene))
    print(on.genes)
  }
  return(on.genes)
}

# Get list of "on" genes
chunckedgenes <- split(intersected.genes, ceiling(seq_along(intersected.genes)/20))
on.progens.18 <- filterAvgGeneExpression(chunckedgenes, daniocell.progens.18)

# filter this list of genes by median expression
# Chunk the genes for functionality
chunked.genes.18 <- split(on.progens.18, ceiling(seq_along(on.progens.18)/20))
on.progens.18.med <- filterMedianGeneExpression(chunked.genes.18, daniocell.progens.18)

on.optic.genes <- filterAvgGeneExpression(chunckedgenes, daniocell.optic)

# # only include genes with mean gene expression > 0.3
# for (var in chunckedgenes) {
#   print(var)
#   optic.cup <- generateData(var, daniocell.optic)
#   
#   mean1 <- optic.cup %>% 
#     group_by(gene, bin) %>% 
#     summarise(mean = mean(geneExpression))
#   
#   mean1 <- mean1 %>% filter(mean > 0.3)
#   print(unique(mean1$gene))
#   on.optic.genes <- c(on.optic.genes, unique(mean1$gene))
#   print(on.optic.genes)
# }

# Export list of on optic genes
write_tsv(tibble(on.optic.genes), "on_optic_genes.txt")
## PRINTING SPECIFIC GENE GRAPHS
#Run til get genes of interest before this code
# counter = 0
# 
# for (Group in chunked.group.genes){
#   print(Group)
#   filtered <- generateData(Group, daniocell.optic)
#   boxplot1 <- ggplot(filtered, aes(x=gene, y=geneExpression, fill=bin))+
#    geom_boxplot()+
#     # geom_point(aes(color=identity.super))+
#    geom_hline(yintercept = 0.3, linetype = "dashed", col = "red") +
#    stat_summary(fun.y = mean, aes(group = bin),position = position_dodge(0.75), geom = "point", color = "lightblue", size=2) +
#    scale_fill_discrete(labels= c("bin.1" = "12-24hpf", "bin.2" = "24-36hpf", "bin.3" = "36-48hpf",
#                                   "bin.4" = "48-60hpf",
#                                   "bin.5" = "60-72hpf", "bin.6" = "72-84hpf",
#                                   "bin.7" = "84-96hpf", "bin.8" = "96-108hpf",
#                                   "bin.9" = ">108hpf"))
#   counter = counter + 1
#   outfile <- paste( Group.num, "-", counter, ".png")
#   ggsave( outfile, boxplot1, device="png", height=1600, width=3200, units = "px")
#   print(outfile)
# }

## LOOPING THROUGH AND SAVING PLOTS

# Make sure to reset counter at 0 when starting anew
saveBoxplot <- function(chunked.on.genes, daniocell.celltype, file_prefix) {
  # Reset counter
  counter = 0
  # Loop through list of "chunked on genes" and save boxplots.
  for (var in chunked.on.genes) {
    print(var)
    filtered <- generateData(var, daniocell.celltype)
    boxplot1 <- ggplot(filtered, aes(x=gene, y=geneExpression, fill=bin))+
      geom_boxplot()+
      # geom_point(aes(color=identity.super))+
      geom_hline(yintercept = 0.3, linetype = "dashed", col = "red") +
      stat_summary(fun.y = mean, aes(group = bin),position = position_dodge(0.75), geom = "point", color = "lightblue", size=2) +
      scale_fill_discrete(labels= c("bin.1" = "12-24hpf", "bin.2" = "24-36hpf", "bin.3" = "36-48hpf",
                                    "bin.4" = "48-60hpf",
                                    "bin.5" = "60-72hpf", "bin.6" = "72-84hpf",
                                    "bin.7" = "84-96hpf", "bin.8" = "96-108hpf",
                                    "bin.9" = ">108hpf"))
    
    counter = counter + 1
    outfile_path <- paste0(paste0(file_prefix, counter), ".png")
    ggsave(outfile_path, boxplot1, device="png", height=1600, width=3200, units = "px")
    print(paste("Saved file", outfile_path))
  }
}

# Generate and save all boxplots for eye.18
chunked.on.genes.18 <- split(on.progens.18.med, ceiling(seq_along(on.progens.18)/8))
saveBoxplot(chunked.on.genes.18, daniocell.progens.18, "onGenes18_")

# Generate and save all boxplots for optic cup
chunked.on.genes.optic <- split(on.optic.genes, ceiling(seq_along(on.optic.genes)/8))
saveBoxplot(chunked.on.genes.optic, daniocell.optic, "onGenesOptic_")

# Generate boxplots for groups for any cluster type

Group.num <- "rpgroup6" #Specify group number
Group.genes <- scan(Group.num, what = "", sep = "\n")

chunked.group.genes <- split(Group.genes, ceiling(seq_along(Group.genes)/8))
saveBoxplot(chunked.group.genes, daniocell.progens.18, paste0(paste0("onGenes18_", Group.num, "_")))
# interchange the tissue type input name and the file output name


# for (var in chunked.on.genes) {
#   print(var)
#   filtered <- generateData(var, daniocell.optic)
#   boxplot1 <- ggplot(filtered, aes(x=gene, y=geneExpression, fill=bin))+
#     geom_boxplot()+
#     # geom_point(aes(color=identity.super))+
#     geom_hline(yintercept = 0.3, linetype = "dashed", col = "red") +
#     stat_summary(fun.y = mean, aes(group = bin),position = position_dodge(0.75), geom = "point", color = "lightblue", size=2) +
#     scale_fill_discrete(labels= c("bin.1" = "12-24hpf", "bin.2" = "24-36hpf", "bin.3" = "36-48hpf",
#                                   "bin.4" = "48-60hpf",
#                                   "bin.5" = "60-72hpf", "bin.6" = "72-84hpf",
#                                   "bin.7" = "84-96hpf", "bin.8" = "96-108hpf",
#                                   "bin.9" = ">108hpf"))
#   
#   counter = counter + 1
#   outfile_path <- paste0(paste0("onGenesBoxplot", counter), ".png")
#   ggsave(outfile_path, boxplot1, device="png", height=1600, width=3200, units = "px")
#   print(paste("Saved file", outfile_path))
# }

## Generate Boxplot for one group of genes
# filtered <- generateData(chunked.on.genes[[1]], daniocell.optic)
# boxplot1 <- ggplot(filtered, aes(x=gene, y=geneExpression, fill=bin))+
#   geom_boxplot()+
#   # geom_point(aes(color=identity.super))+
#   geom_hline(yintercept = 0.3, linetype = "dashed", col = "red") +
#   stat_summary(fun.y = mean, aes(group = bin),position = position_dodge(0.75), geom = "point", color = "lightblue", size=2) +
#   scale_fill_discrete(labels= c("bin.1" = "12-24hpf", "bin.2" = "24-36hpf", "bin.3" = "36-48hpf",
#                                 "bin.4" = "48-60hpf",
#                                 "bin.5" = "60-72hpf", "bin.6" = "72-84hpf",
#                                 "bin.7" = "84-96hpf", "bin.8" = "96-108hpf",
#                                 "bin.9" = ">108hpf"))
# 
# boxplot1
# counter <- (counter+1)
# outfile_path <- paste0(paste0("onGenesBoxplot", counter), ".png")
# ggsave(outfile_path, boxplot1, device="png", height=1600, width=3200, units = "px")


#bmp7btest <- generateData("bmp7b")
#bmp7b <- inner_join(bmp7btest, eye.clust.names, by = "cell")
#retinal.progen <- bmp7b %>% filter(identity.super == "retinal progenitors")
#optic.cup <- bmp7b %>% filter(identity.super == "optic cup")
# perCells <- paste0(round(count(Agg0)/count(Agg)*100, 2), '% of eye cells')

## FILTERING BY AVERAGE GENE EXPRESSION-----------

mean1 <- optic.cup %>% 
  group_by(gene, bin) %>% 
  summarise(mean = mean(geneExpression))

mean1 <- mean1 %>% filter(mean > 0.3)
unique(mean1$gene)

filtered <- optic.cup %>% filter(gene %in% unique(mean1$gene))

## PLOTTING BOXPLOTS---------------------------------
ggplot(filtered, aes(x=gene, y=geneExpression, fill=bin))+
  geom_boxplot()+
  # geom_point(aes(color=identity.super))+
  geom_hline(yintercept = 0.3, linetype = "dashed", col = "red") +
  stat_summary(fun.y = mean, aes(group = bin),position = position_dodge(0.75), geom = "point", color = "lightblue", size=2) +
  scale_fill_discrete(labels= c("bin.1" = "12-24hpf", "bin.2" = "24-36hpf", "bin.3" = "36-48hpf",
                                "bin.4" = "48-60hpf",
                                "bin.5" = "60-72hpf", "bin.6" = "72-84hpf",
                                "bin.7" = "84-96hpf", "bin.8" = "96-108hpf",
                                "bin.9" = ">108hpf"))
  # scale_color_discrete(color=bmp7b$identity.super)
  # facet_wrap(~gene)
  # annotate("text", x=1.3, y = 4, label = perCells)

# intersecting gene lists----------------------
partofrp <- read_csv("C:/Users/ayahs/Downloads/partofrp", col_names=FALSE)
partofoc <- read_csv("C:/Users/ayahs/Downloads/partofoc", col_names=FALSE)
partOfBoth <- intersect(partofrp, partofoc)

## PLOTTING GENE EXPRESSION ------------

# Plot a gene on Daniocell
FeaturePlot(daniocell, "rx3")

# Plot a gene on Daniocell, eye
FeaturePlot(daniocell.eye, "rx3")

# Make a violin plot to show the distribution of gene expression across several clusters
VlnPlot(daniocell.eye, "olig2")
VlnPlot(daniocell.eye, c("rx3", "olig2"))
VlnPlot(daniocell.eye, c("rx3", "olig2"), idents = c("optic field", "optic cup", "retinal progenitors", "retinal ganglion cells"))

# Make a dot plot to show expression summarized for many genes across clusters
DotPlot(daniocell.eye, c("rx3", "olig2", "sox2", "atoh7", "col15a1b", "vsx2"))
DotPlot(daniocell.eye, c("rx3", "olig2", "sox2", "atoh7", "col15a1b", "vsx2"), scale = F)

# Get number of observations per bin
observationsPerBin <- function(daniocell.celltype) {
  after.12 <- WhichCells(daniocell.celltype, expression = hpf >= 12)
  after.14 <- WhichCells(daniocell.celltype, expression = hpf >= 14)
  after.24 <- WhichCells(daniocell.celltype, expression = hpf >= 24)
  after.36 <- WhichCells(daniocell.celltype, expression = hpf >= 36)
  after.48 <- WhichCells(daniocell.celltype, expression = hpf >= 48)
  after.60 <- WhichCells(daniocell.celltype, expression = hpf >= 60)
  after.72 <- WhichCells(daniocell.celltype, expression = hpf >= 72)
  after.84 <- WhichCells(daniocell.celltype, expression = hpf >= 84)
  after.96 <- WhichCells(daniocell.celltype, expression = hpf >= 96)
  bin.9 <- WhichCells(daniocell.celltype, expression = hpf >= 108)
  bin.8 <- setdiff(after.96, bin.9)
  bin.7 <- setdiff(after.84, after.96)
  bin.6 <- setdiff(after.72, after.84)
  bin.5 <- setdiff(after.60, after.72)
  bin.4 <- setdiff(after.48, after.60)
  bin.3 <- setdiff(after.36, after.48)
  bin.2 <- setdiff(after.24, after.36)
  bin.1 <- setdiff(after.12, after.24)

  print(paste("Bin 1", length(bin.1)))
  print(paste("Bin 2", length(bin.2)))
  print(paste("Bin 3", length(bin.3)))
  print(paste("Bin 4", length(bin.4)))
  print(paste("Bin 5", length(bin.5)))
  print(paste("Bin 6", length(bin.6)))
  print(paste("Bin 7", length(bin.7)))
  print(paste("Bin 8", length(bin.8)))
  print(paste("Bin 9", length(bin.9)))
}

observationsPerBin(daniocell.optic)
observationsPerBin(daniocell.progens.18)
