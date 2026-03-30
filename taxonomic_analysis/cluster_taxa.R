### Load libraries 
library(phyloseq)
library(pheatmap)

### Read in files output by generate_phyloseq.py
labeled.metadata <- read.csv("finished_scripts/labMeta.csv", header=T, row.names = 1, na.strings = c("NA", "ND"))
labeled.metadata$Name <- rownames(labeled.metadata)
labeled.SampleData = sample_data(labeled.metadata)

labeled.taxonomy.table <- read.csv("finished_scripts/labTax.csv", header=T, row.names = 1, na.strings = c("NA", "ND"))

labeled.OTU.table <- read.csv("finished_scripts/labOTU.csv", row.names = 1, header=T, check.names=FALSE)
AE.OTU.table <- read.csv("finished_scripts/totAE.csv", row.names = 1, header=T, check.names=FALSE)

### Values = labeled spectral counts
labeled.taxonomy.mat <- tax_table(as.matrix(labeled.taxonomy.table))
labeled.OTU.mat <- otu_table(as.matrix(labeled.OTU.table), taxa_are_rows = TRUE)

### Values = average enrichment 
AE.OTU.mat <- otu_table(as.matrix(AE.OTU.table), taxa_are_rows = TRUE)

### Generate phyloseq objects
labeledSCs.phyloseq <- phyloseq(labeled.taxonomy.mat, labeled.OTU.mat, labeled.SampleData)
totalAE.phyloseq <- phyloseq(labeled.taxonomy.mat, AE.OTU.mat, labeled.SampleData)

### Subset phyloseq objects so only cecum samples are included
labeledSC.phyloseq.cecum <- subset_samples(labeledSCs.phyloseq, SampleType == "C")
totalAE.phyloseq.cecum <- subset_samples(totalAE.phyloseq, SampleType == "C")

### Read in CSV containing signifincantly labeled taxa
sigTaxa <- read.csv('Significantly_Labeled_Taxa.csv', row.names = 1)
sigTaxa <- row.names(sigTaxa)

### Cluster based on labeled spectral counts
labeledSC.phyloseq.taxmat <- as(otu_table(labeledSC.phyloseq.cecum), "matrix")
labeledSC.phyloseq.SIG.taxmat <- labeledSC.phyloseq.taxmat[rownames(labeledSC.phyloseq.taxmat) %in% sigTaxa, ]

labeledSC.phyloseq.taxmat.dist <- dist(labeledSC.phyloseq.SIG.taxmat, method = "euclidean")
labeledSC.phyloseq.taxmat.clust <- hclust(labeledSC.phyloseq.taxmat.dist, method = "ward.D2")

rownames(labeledSC.phyloseq.SIG.taxmat) <- gsub("^[a-z]__", "", rownames(labeledSC.phyloseq.SIG.taxmat))
time <- as.numeric(gsub(".*C(\\d+).*", "\\1", colnames(labeledSC.phyloseq.SIG.taxmat)))
labeledSC.phyloseq.SIG.taxmat <- labeledSC.phyloseq.SIG.taxmat[, order(time), drop = FALSE]

pheatmap(
  labeledSC.phyloseq.SIG.taxmat,
  color = hcl.colors(50, "YlOrRd"),
  cluster_rows = labeledSC.phyloseq.taxmat.clust,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  border_color = "white",
  main = 'Labeled Spectral Count'
)

### Cluster based on average enrichment of labeled PSMs
totalAE.phyloseq.taxmat <- as(otu_table(totalAE.phyloseq.cecum), "matrix")
totalAE.phyloseq.SIG.taxmat <- totalAE.phyloseq.taxmat[rownames(totalAE.phyloseq.taxmat) %in% sigTaxa, ]

totalAE.phyloseq.taxmat.dist <- dist(totalAE.phyloseq.SIG.taxmat, method = "euclidean")
totalAE.phyloseq.taxmat.clust <- hclust(totalAE.phyloseq.taxmat.dist, method = "ward.D2")

rownames(totalAE.phyloseq.SIG.taxmat) <- gsub("^[a-z]__", "", rownames(totalAE.phyloseq.SIG.taxmat))
time <- as.numeric(gsub(".*C(\\d+).*", "\\1", colnames(totalAE.phyloseq.SIG.taxmat)))
totalAE.phyloseq.SIG.taxmat <- totalAE.phyloseq.SIG.taxmat[, order(time), drop = FALSE]

pheatmap(
  totalAE.phyloseq.SIG.taxmat,
  color = hcl.colors(50, "YlOrRd"),
  cluster_rows = totalAE.phyloseq.taxmat.clust,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  fontsize_row = 10,
  border_color = "white",
  main = 'Average Enrichment'
)

