### Load libraries 
library(phyloseq)
library(vegan)
library(ggplot2)
library(compositions)

#### COMPARE LAB VS. UNLABE
unlabTAXONMetadata <- read.csv("metadata_table_unlab.csv", header=T, row.names = 1, na.strings = c("NA", "ND"))
labTAXONMetadata <- read.csv("metadata_table_lab.csv", header=T, row.names = 1, na.strings = c("NA", "ND"))
unlabTAXONMetadata$Name <- rownames(unlabTAXONMetadata)
labTAXONMetadata$Name <- rownames(labTAXONMetadata)

unlabeledTAXONTaxonomyTable <- read.csv("unlabeled_taxonomy_table.csv", header=T, row.names = 1, na.strings = c("NA", "ND"))
labeledTAXONTaxonomyTable <- read.csv("labeled_taxonomy_table.csv", header=T, row.names = 1, na.strings = c("NA", "ND"))
unlabeledTAXONOTUTable <- read.csv("unlabeled_otu_table.csv", row.names = 1, header=T, check.names=FALSE)
labeledTAXONOTUTable <- read.csv("labeled_otu_table.csv", row.names = 1, header=T, check.names=FALSE)


labeledTAXONTaxonomyMatrix <- tax_table(as.matrix(labeledTAXONTaxonomyTable))
labeledTAXONOTUMatrix <- otu_table(as.matrix(labeledTAXONOTUTable), taxa_are_rows = TRUE)

unlabeledTAXONTaxonomyMatrix <- tax_table(as.matrix(unlabeledTAXONTaxonomyTable))
unlabeledTAXONOTUMatrix <- otu_table(as.matrix(unlabeledTAXONOTUTable), taxa_are_rows = TRUE)

unlabeledTAXONSampleData = sample_data(unlabTAXONMetadata)
labeledTAXONSampleData = sample_data(labTAXONMetadata)

labeledPhyloseqObj.taxon <- phyloseq(labeledTAXONTaxonomyMatrix, labeledTAXONOTUMatrix, labeledTAXONSampleData)
unlabeledPhyloseqObj.taxon <- phyloseq(unlabeledTAXONTaxonomyMatrix, unlabeledTAXONOTUMatrix, unlabeledTAXONSampleData)
labeledPhyloseqObj.taxon.cecum <- subset_samples(labeledPhyloseqObj.taxon, Type == "Cecum")
unlabeledPhyloseqObj.taxon.cecum <- subset_samples(unlabeledPhyloseqObj.taxon, Type == "Cecum")

### NMDS ############
labeledNMDS.taxon.cecum <- ordinate(labeledPhyloseqObj.taxon.cecum, "NMDS", "bray")
plot_ordination(labeledPhyloseqObj.taxon.cecum, labeledNMDS.taxon.cecum, color="Group", title = 'NMDS_labeled_cecum') +
  geom_text(aes(label = Name), size = 3)

unlabeledNMDS.taxon.cecum <- ordinate(unlabeledPhyloseqObj.taxon.cecum, "NMDS", "bray")
plot_ordination(unlabeledPhyloseqObj.taxon.cecum, unlabeledNMDS.taxon.cecum, color="Group", title = 'NMDS_unlabeled_cecum') +
  geom_text(aes(label = Name), size = 3, vjust = -0.8)

##### PERMANOVA SIGNIFICANCE TESTING ON CLR-TRANSFORMED DATA
labeledOTUs.taxon.pseudocounts <- labeledTAXONOTUMatrix + 1
labeledOTUs.taxon.clrTransformed <- t(clr(acomp(t(labeledOTUs.taxon.pseudocounts))))
labeledOTUTable.taxon.clrTransformed <- otu_table(labeledOTUs.taxon.clrTransformed, taxa_are_rows = TRUE)
labeledPhyloseqObj.taxon.clr <- phyloseq(labeledTAXONTaxonomyMatrix, labeledOTUTable.taxon.clrTransformed, labeledTAXONSampleData)
labeledPhyloseqObj.taxon.clr.stool <- subset_samples(labeledPhyloseqObj.taxon.clr, Type == "Stool")
labeledPhyloseqObj.taxon.clr.cecum <- subset_samples(labeledPhyloseqObj.taxon.clr, Type == "Cecum")

### LABELED
labeledPERMANOVAMetadata.taxon.clr.cecum <- as.data.frame(labeledPhyloseqObj.taxon.clr.cecum@sam_data@.Data) #create a data.frame with the metadata from the phyloseq object; I don't use the "meta" object created above, because sometimes the sample counts aren't the same
rownames(labeledPERMANOVAMetadata.taxon.clr.cecum) <- labeledPhyloseqObj.taxon.clr.cecum@sam_data@row.names #label the rownames
colnames(labeledPERMANOVAMetadata.taxon.clr.cecum) <- labeledPhyloseqObj.taxon.clr.cecum@sam_data@names #label the colnames

labeledPERMANOVAMetadata.taxon.clr.O.cecum <- labeledPERMANOVAMetadata.taxon.clr.cecum[order(rownames(labeledPERMANOVAMetadata.taxon.clr.cecum)),] #adonis2 needs the metadata table to be ordered the same as the distance matrix
labeledBrayDist.taxon.clr.cecum <- vegdist(t(labeledPhyloseqObj.taxon.clr.cecum@otu_table), method = "euclidean", diag = TRUE) #create a distance matrix from the asv count table of the phyloseq object
labeledBrayDistMatrix.taxon.clr.cecum <- as.matrix(labeledBrayDist.taxon.clr.cecum)
labeledPERMANOVAMetadata.taxon.clr.O.cecum <- labeledPERMANOVAMetadata.taxon.clr.cecum[labels(labeledBrayDist.taxon.clr.cecum), , drop = FALSE]

identical(row.names(labeledBrayDistMatrix.taxon.clr.cecum), row.names(labeledPERMANOVAMetadata.taxon.clr.cecum)) #check that rownames match between the two objects
identical(rownames(labeledPERMANOVAMetadata.taxon.clr.O.cecum), labels(labeledBrayDist.taxon.clr.cecum))

adonis2(labeledBrayDist.taxon.clr.cecum ~ Time, data = labeledPERMANOVAMetadata.taxon.clr.O.cecum, permutations = 999) 

### UNLABELED
unlabeledOTUs.taxon.pseudocounts <- unlabeledTAXONOTUMatrix + 1
unlabeledOTUs.taxon.clrTransformed <- t(clr(acomp(t(unlabeledOTUs.taxon.pseudocounts))))
unlabeledOTUTable.taxon.clrTransformed <- otu_table(unlabeledOTUs.taxon.clrTransformed, taxa_are_rows = TRUE)
unlabeledPhyloseqObj.taxon.clr <- phyloseq(unlabeledTAXONTaxonomyMatrix, unlabeledOTUTable.taxon.clrTransformed, unlabeledTAXONSampleData)
unlabeledPhyloseqObj.taxon.clr.stool <- subset_samples(unlabeledPhyloseqObj.taxon.clr, Type == "Stool")
unlabeledPhyloseqObj.taxon.clr.cecum <- subset_samples(unlabeledPhyloseqObj.taxon.clr, Type == "Cecum")

unlabeledPERMANOVAMetadata.taxon.clr.cecum <- as.data.frame(unlabeledPhyloseqObj.taxon.clr.cecum@sam_data@.Data) #create a data.frame with the metadata from the phyloseq object; I don't use the "meta" object created above, because sometimes the sample counts aren't the same
rownames(unlabeledPERMANOVAMetadata.taxon.clr.cecum) <- unlabeledPhyloseqObj.taxon.clr.cecum@sam_data@row.names #label the rownames
colnames(unlabeledPERMANOVAMetadata.taxon.clr.cecum) <- unlabeledPhyloseqObj.taxon.clr.cecum@sam_data@names #label the colnames

unlabeledPERMANOVAMetadata.taxon.clr.O.cecum  <- unlabeledPERMANOVAMetadata.taxon.clr.cecum[order(rownames(unlabeledPERMANOVAMetadata.taxon.clr.cecum)),] #adonis2 needs the metadata table to be ordered the same as the distance matrix
unlabeledBrayDist.taxon.clr.cecum <- vegdist(t(unlabeledPhyloseqObj.taxon.clr.cecum@otu_table), method = "euclidean", diag = TRUE) #create a distance matrix from the asv count table of the phyloseq object
unlabeledBrayDistMatrix.taxon.clr.cecum <- as.matrix(unlabeledBrayDist.taxon.clr.cecum)
unlabeledPERMANOVAMetadata.taxon.clr.O.cecum <- unlabeledPERMANOVAMetadata.taxon.clr.cecum[labels(unlabeledBrayDist.taxon.clr.cecum), , drop = FALSE]

identical(row.names(unlabeledBrayDistMatrix.taxon.clr.cecum), row.names(unlabeledPERMANOVAMetadata.taxon.clr.O.cecum)) #check that rownames match between the two objects
identical(rownames(unlabeledPERMANOVAMetadata.taxon.clr.O.cecum), labels(unlabeledBrayDist.taxon.clr.cecum))

adonis2(unlabeledBrayDist.taxon.clr.cecum ~ Time, data = unlabeledPERMANOVAMetadata.taxon.clr.O.cecum, permutations = 999)
