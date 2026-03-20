#!/usr/bin/env python3

'''
taxonomy_bar_chart.py

Purpose: 
    Generate bar chart displaying proportional labeled and unlabeled spectral counts of most abundant genera

Inputs:
    - Labeled proteome output by Percolator (TSV)
    - Unlabeled proteome parsed from Sipros output (TSV)
    - File names --> sample names lookup table (CSV)
    - MGYG proteome metadata with taxon IDs --> taxon lineage 

Outputs:
    Taxonomy bar chart

Usage:
    python taxonomy_bar_chart.py \
        -l [labeled proteome] \
        -u [unlabeled proteome] \
        -n [sample names lookup table] \
        -m [MGYG proteome metadata]
        
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

class parseMetadata():
    """
    Parse sample naming metadata and MGYG metadata 

    Generate dictionaries to query sample names, label status, genera names, 
    and entire lineage of each genus

    Attributes:
        namesDict (file path) : path to sample naming metadata 
        MGYGMeta (file path) : path to MGYG taxonomy metadata 
    """

    def sampleMetadata(self, namesDictIn):
        """
        Generate lookup dictionaries containing metadata for all samples

        Parameters
        ----------
        namesDictIn : file path
            CSV file contains lookup table to convert MS/MS file name --> experimental sample name 

        Returns 
        -------
        groupDict : dict
            Dictionary to convert MS/MS file name to treatment group ID
        
        sampStatDict : dict 
            Dictionary to lookup if sample is unlabeled or labeled 

        Notes
        -----
        Samples from this experiment are named like [Sample type][Time point].[Replicate]    
        """
        sampleLookup = pd.read_csv(namesDictIn)
        sampleLookupDict = sampleLookup.to_dict(orient = 'index')
        sampleDict = {}
        groupDict = {}
        sampStatDict = {}
        orderDict = {}
        for sampleID, sampleName in sampleLookupDict.items():
            group = sampleName['SampleName'].split('.')
            sampleDict[sampleName['FileName']] = sampleName['SampleName']
            groupDict[sampleName['SampleName']] = group[0]
            orderDict[sampleName['SampleName']] = int(group[1])
            if '0' in sampleName['SampleName'].split('.')[0]:
                sampStatDict[sampleName['SampleName']] = 'Unlabeled'
            else:
                sampStatDict[sampleName['SampleName']] = 'Labeled'
        return sampleDict, groupDict, sampStatDict, orderDict
    

    def parseMGYGData(self, metadata):
        """
        Parse MGYG metadata file into dictionary that can be used to retrieve genus names for each taxon ID

        Parameters
        ----------
        metadata : file path to MGYG metadata CSV
            CSV file contains many fields, here we only read in the taxon ID and full lineage with taxon names

        Returns 
        -------
        lineageDict : dict 
            Lookup dict to get genus name for each taxon ID
        """
        metadataDf = pd.read_csv(metadata, sep = '\t', header = 0, usecols = ['Genome', 'Lineage'])
        lineageDict = {}
        for isolate, lineage in metadataDf.itertuples(index = False):
            splitLineage = lineage.split(';')
            ### splitLineage is a list that can be indexed to get any taxonomic rank in the lineage
            ### Index 5 corresponds to genera 
            genus = splitLineage[5].split('__')[1]
            lineageDict[isolate] = genus
        return lineageDict

class computeTaxonSpectralCounts():
    """
    Parse data and plot proportional spectral counts of most abundant taxa

    This class parses the labeled and unlabeled proteome files, calculates the proportional spectral counts
    of most abundant genera, and then plots those values as a bar chart
    """

    def __init__(self, sDict, gDict, stDict, oDict, linDict):
        self.sDict = sDict
        self.gDict = gDict
        self.stDict = stDict
        self.oDict = oDict
        self.linDict = linDict

    def parseAbundance(self, df, status):
        """
        Parse labeled or unlabeled proteome data

        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame containing PSM data

        status : list
            String to indicate which types of samples (Labeled or Unlabeled) 
            should be included in the resulting plot 

        Returns
        -------
        transformSCDf : pandas.DataFrame
            DataFrame with Protein ID as index, spectral counts as values, and 
            hierarchical columns with information about time point and sample 

        Notes
        -----
        All samples are included in unlabeled proteome plot because there are 
        unlabeled PSMs in all samples. Only labeled samples are included in the labeled
        proteome plot because any labeled detections in unlabeled samples are known
        to be false positives
        """
        abundData = []
        for psm, protein, sample in df.itertuples(index = False):
            stripProtein = protein.lstrip('{').rstrip('}')
            sampleName = self.sDict.get(sample)
            sampleStatus = self.stDict.get(sampleName)
            ### Only plot cecum samples
            if 'C' in sampleName and sampleStatus in status:
                if stripProtein.startswith('MGYG'):
                    splitProtein = stripProtein.split(',')
                    ### Append first protein in list of proteins 
                    ### if it's a degenerate protein, first entry is saved
                    abundData.append([splitProtein[0], sampleName])
    
        abundDataDf = pd.DataFrame(abundData).rename(columns = {0: 'Protein', 1: 'Sample'})

        ### Make sure samples are in chronological order
        abundDataDf['Order'] = abundDataDf['Sample'].map(self.oDict)
        abundDataDf = abundDataDf.sort_values(by = 'Order').drop('Order', axis = 1)

        ### Use groupby function to compute spectral count of proteins, assume each PSM represents one spectral count
        gbProtein = abundDataDf.groupby(['Protein', 'Sample'], sort = False).size().reset_index()
        gbProtein['Group'] = gbProtein['Sample'].map(self.gDict)
        transformSCDf = pd.pivot_table(gbProtein, index = 'Protein', columns = ['Sample', 'Group'], aggfunc = 'sum', sort = False).fillna(0)
        transformSCDf.columns = transformSCDf.columns.droplevel()
        return transformSCDf
    
    def insertTaxonomy(self, scDf):
        """
        Insert genera names associated with each protein into dataframe

        Parameters
        ----------
        scDf : pandas.DataFrame
            DataFrame with samples as columns, index as protein IDs, and values as spectral counts

        Returns
        -----
        gbTaxon : pandas.DataFrame
            Spectral counts summarized at the genus level
        """
        taxonomyData = {}
        for proteinID, sampleData in scDf.iterrows():
            taxonID = proteinID.split('_')[0]
            taxon = self.linDict.get(taxonID)
            taxonomyData[proteinID] = taxon
        scDf = scDf.reset_index()
        scDf['Taxon'] = scDf['Protein'].map(taxonomyData)
        scDf = scDf.set_index('Taxon').drop('Protein', axis = 1, level = 0).reset_index()
        gbTaxon = scDf.groupby('Taxon').sum()
        return gbTaxon

    def computeProportions(self, spectralCountDf):
        """
        Normalize raw spectral count data by transforming them into proportions

        Parameters
        ----------
        spectralCountDf : pandas.DataFrame
            Raw spectral counts summarized at the genus level

        Returns
        -------
        proportionsDf : pandas.DataFrame
            DataFrame containing proportional spectral counts of most abundant genera in labeled 
            and unlabeled proteomes
        """
        totalAbundances = spectralCountDf.sum(axis = 1)
        abundantTaxa = totalAbundances.nlargest(8).index
        mostAbundantTaxaGb = spectralCountDf.loc[abundantTaxa]

        ### Summarize proportion of spectral counts assigned to less abundant taxa as "Other"
        others = spectralCountDf.loc[~spectralCountDf.index.isin(abundantTaxa)].sum()
        otherRowData = others.to_dict()
        mostAbundantTaxaGb.loc['Other'] = otherRowData
        proportionsDf = mostAbundantTaxaGb/mostAbundantTaxaGb.sum()
        return proportionsDf
    
    def assignColors(self, df1, df2):
        """
        Generate custom colormap by assigning bar chart colors to each taxon that
        will be plotted

        Parameters
        ----------
        df1 : pandas.DataFrame
            DataFrame containing labeled proteome data

        df2 : pandas.DataFrame
            DataFrame containing unlabeled proteome data

        Returns
        -------
        colorMap : dict
            Dictionary to be used to map genus name to corresponding color on plot
        """
        allTaxa = df1.index.values.tolist() + df2.index.values.tolist()
        uniqueTaxa = pd.Series(allTaxa).drop_duplicates().values.tolist()
        colors = [
                "lightblue",
                "steelblue",
                "yellowgreen",
                "forestgreen",
                "navy",
                "firebrick",
                "navajowhite",
                "darkorange",
                "lightgrey",
                "slategray",
                "darkmagenta",
                "rosybrown",
                "teal",
                "tan",
                "crimson",
            ]
        colorMap = {}
        for taxon, c in zip(uniqueTaxa, colors):
            colorMap[taxon] = c
        return colorMap

    def plotTaxa(self, plotDf, cmap, axis, plotTitle):
        """
        Plot proportional spectral counts of genera in labeled and unlabeled proteomes

        Parameters
        ----------
        plotDf : pandas.DataFrame
            DataFrame with data from labeled or unlabeled proteome
        
        cmap : dict
            Dictionary to be used to map genus name to corresponding color on plot
        
        axis : matplotlib.axes._axes.Axes
            Axis to use to plot data
        
        plotTitle : str
            Subplot title 
        """
        gap = 2
        xPos = 0
        xTicks = []
        sampleLabs = []
        for i, (tGroup, tGroupData) in enumerate(plotDf.T.groupby(level="Group", sort = False)):
            for sample, row in tGroupData.iterrows():
                bottoms = 0
                sampleLabs.append(sample[0])
                for plotTax, taxonAbundance in row.reset_index().itertuples(index = False):
                    barColor = cmap.get(plotTax)
                    axis.bar(xPos, taxonAbundance, bottom = bottoms, width = 0.9, color = barColor, label = plotTax)
                    bottoms += taxonAbundance
                xTicks.append(xPos)
                xPos += 1
            xPos += gap
        axis.set_xticks(xTicks)
        axis.set_xticklabels(sampleLabs, rotation = 90, fontsize = 11)
        axis.set_ylabel("Proportion")
        axis.set_xlabel("Sample")
        axis.set_title(plotTitle, pad=14)
            
    def plotLegend(self, axes):
        """
        Plot legend 

        Parameters
        ----------
        axes : np.array
            Array containing both axes that will be plotted
        """
        handles, labels = [], []
        for a in axes.ravel():
            h, l = a.get_legend_handles_labels()
            for hand, lab in zip(h, l):
                if lab and not lab.startswith('_') and lab not in labels:
                    handles.append(hand)
                    labels.append(lab)
        
        ### Make it so "Other" is last in legend
        if "Other" in labels:
            otherIDX = labels.index("Other")
            oh = handles.pop(otherIDX)
            ol = labels.pop(otherIDX)
            handles.append(oh)
            labels.append(ol)

        axes[1].legend(handles, labels, bbox_to_anchor=(1.01,0.98), title='Genus', fontsize = 10, frameon = False)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--labeledProteome')
    parser.add_argument('-u', '--unlabeledProteome')
    parser.add_argument('-n', '--namesDict')
    parser.add_argument('-m', '--metadata')
    args = parser.parse_args()

    labeledData = pd.read_csv(args.labeledProteome, sep = '\t', header = 0, usecols = [0, 25, 26])
    unlabeledData = pd.read_csv(args.unlabeledProteome, sep = '\t', header = 0, usecols = [0, 26, 27])

    parseMeta = parseMetadata()
    sampleDict, groupDict, statusDict, ordDict = parseMeta.sampleMetadata(args.namesDict)
    taxonomyLookupDict = parseMeta.parseMGYGData(args.metadata)

    parseSCs = computeTaxonSpectralCounts(sampleDict, groupDict, statusDict, ordDict, taxonomyLookupDict)
    labSpectralCountDf = parseSCs.parseAbundance(labeledData, ['Labeled'])
    unlabSpecralCountDf = parseSCs.parseAbundance(unlabeledData, ['Labeled', 'Unlabeled'])
    
    labTaxaSCDf = parseSCs.insertTaxonomy(labSpectralCountDf)
    unlabTaxaSCDf = parseSCs.insertTaxonomy(unlabSpecralCountDf)

    labeledTaxonProportions = parseSCs.computeProportions(labTaxaSCDf)
    unlabeledTaxonProportions = parseSCs.computeProportions(unlabTaxaSCDf)

    colorMapDict = parseSCs.assignColors(labeledTaxonProportions, unlabeledTaxonProportions)
    
    fig, ax = plt.subplots(1, 2, figsize=(15, 7))
    fig.subplots_adjust(left=0.069, top = 0.9, right = 0.8, bottom = 0.22, wspace = 0.3)

    parseSCs.plotTaxa(unlabeledTaxonProportions, colorMapDict, ax[0], "Unlabeled Spectral Counts")
    parseSCs.plotTaxa(labeledTaxonProportions, colorMapDict, ax[1], "Labeled Spectral Counts")

    parseSCs.plotLegend(ax)
    plt.show()

if __name__ == "__main__":
    main()
