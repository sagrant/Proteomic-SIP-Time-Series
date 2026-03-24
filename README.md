# Proteomic SIP Replicated Time Series Experiment 
This repository contains code and metadata required to visualize and analyze data produced by a Proteomic SIP experiment examining <sup>13</sup>C algal protein assimilation in the mouse gut microbiome.

<p align="center">
    <img src="images/algal_protein_experiment_overview.svg" width="700">

The gut microbiome is closely related to host health and plays a major role in host nutrition by assisting in the metabolism of dietary substrates that the host ingests. The ecological processes that dictate how dietary nutrients affect the gut microbiome’s ecology, and furthermore, host health, are poorly understood. Proteomic Stable Isotope Probing (SIP) is a method that can link an isotopically-labeled substrate to the organisms that assimilate it, thereby identifying the organisms and putative ecological processes involved in the substrate’s degradation. Proteins synthesized by microbes that assimilate the substrate into anabolic pathways will be enriched with the stable isotope, these are referred to as labeled proteins. Matching labeled mass spectra to known peptide sequences encoded by specific microbial populations enables the identification of active taxa and allows for quantification of their activity and specificity. We applied Proteomic SIP in the mouse gut microbiome to show that this method can be adapted to draw causal inferences about the effect of dietary substrates on the gut microbiome *in vivo*. We showed that Proteomic SIP can be successfully implemented in the mouse gut microbiome. We also demonstrated that this method can detect the biological response of the microbial community to dietary substrates.  

## Objectives
* Count the number of labeled PSMs, peptides, and proteins detected
* Calculate False Positive Rate (FPR)
* Visualize how the proportion of labeled PSMs changes over time and determine if changes are significant
* Visualize distribution of <sup>13</sup>C enrichment over time
* Compare empirical counts of labeled PSMs for detected genera with what would be expected under null conditions
* Identify organisms that assimilated the labeled substrate, and those that did not
* Visualize the relationship between average <sup>13</sup>C enrichment and spectral count for significantly labeled genera

## Requirements and dependencies
* Python 3.13.5
* Pandas 3.0.0
* Numpy 2.4.2
* Matplotlib 3.10.8
* R version 4.5.1
* phyloseq 1.52.0
* vegan 2.7-2
* compositions 2.0-9
* Hotelling 1.0-8

## Directory structure
```
├── images
├── community_analysis
│   ├── count_detected_PSMsPeptidesProteins.py  # Count labeled PSMs, peptides, and proteins
│   ├── calculate_FPR.py    # Calculate dataset FPR
│   ├── proportion_labeled_PSMs.py  # Calculate the proportion of labeled PSMs out of all PSMs 
│   ├── enrichment_distributions.py # Visualize 13C enrichment distributions of all labeled PSMs 
│   ├── generate_phyloseq.py    # Parse data into tables that can be used to make a phyloseq object in R
│   └── NMDS_and_PERMANOVA.R    # Generate NMDS ordination and run PERMANOVA tests
├── taxonomic_analysis
│   ├── null_distributions.py   # Determine which organisms were significantly labeled using null distributions 
│   ├── subset_unlabeled.py # Subset unlabeled proteome from total proteome to be used to generate taxonomy bar chart
│   ├── taxonomy_bar_chart.py   # Visualize how proportional spectral counts of abundant genera change over time 
│   ├── averageEnrichment_spectralCount_bubblePlot.py   # Visualize relationship between average enrichment and spectral count
│   └── cluster_taxa.R  # Cluster genera based on average 13C enrichment or labeled spectral counts
└── phylogenetic and functional analysis
    ├── trait_data.py   # Generate metadata table to be used to make heatmap that accompanies phylogenetic tree
    └── phylogenetic_functions_tree.R   # Generate phylogenetic tree
```