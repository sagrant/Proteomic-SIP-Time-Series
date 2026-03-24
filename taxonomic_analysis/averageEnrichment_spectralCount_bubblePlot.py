#!/usr/bin/env python3

"""
averageEnrichment_spectralCount_visualization.py

Purpose: 
    Visualize the relationship between spectral count and average enrichment for significantly labeled genera.
    This visualization can be used to distinguish between putative resource specialists and generalists.

Inputs:
    - Percolator output (TSV) with all SIP labeled PSMs
    - File names --> sample names lookup table (CSV)
    - MGYG proteome metadata with taxon IDs --> taxon lineage 
    - CSV containing significantly labeled taxa as defined by null distributions

Outputs: 
    - Average enrichment x Spectral count bubble plot

Usage: 
    python averageEnrichment_spectralCount_visualization.py \
        -i [percolator output] \
        -n [sample names lookup table] \
        -m [MGYG metadata] \
        -t [significantly labeled taxa]

Notes:
    [insert justification for thresholds]
    
"""

import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize, LinearSegmentedColormap
from matplotlib.cm import ScalarMappable


class parseSIPData():
    """
    Parses Percolator output into dataframe that can be plotted

    This class integrates sample name metadata, genus name metadata, and SIP-labeled PSMs
    into a datastructure that can be represented as a bubble plot.

    Attributes: 
        SIPdf (pandas.DataFrame): percolator output file data
    """

    def __init__(self, SIPdf):
        self.SIPdf = SIPdf

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
        groupDict = {}
        sampStatDict = {}
        for sampleID, sampleName in sampleLookupDict.items():
            group = sampleName['SampleName'].split('.')
            groupDict[sampleName['FileName']] = group[0]
            if '0' in sampleName['SampleName'].split('.')[0]:
                sampStatDict[group[0]] = 'Unlabeled'
            else:
                sampStatDict[group[0]] = 'Labeled'
        return groupDict, sampStatDict

    def parseInFile(self, treatmentGroupDict, statDict):
        """
        Parse Percolator output into 2 DataFrames with spectral count data or average enrichment data

        Parameters
        ----------
        SIPdf : pandas.DataFrame
            DataFrame containing labeled PSMs, including their enrichment values

        treatmentGroupDict : dict
            Dictionary used to retrieve treatment group based on MS/MS file name

        statDict : dict
            Dictionary used to retrieve status (labeled/unlabeled) of sample based on treatment group

        Returns
        -------
        t_avgEnrichmentDf : pandas.DataFrame
            DataFrame containing average enrichment values, summarized at the protein level
            Columns are samples, index contains protein IDs

        t_spectralCountDf : pandas.DataFrame
            DataFrame containing spectral counts, summarized at the protein level
            Columns are samples, index contains protein IDs
        
        Notes
        -----
        Including degenerate proteins - uses first degenerate entry in list of possible protein matches
        Only includes microbial PSMs
        Excludes all unlabeled samples    
        """
        enrichData = []
        countData = []
        for psm, enrich, ms2, protein, sample in self.SIPdf.itertuples(index = False):
            stripProtein = protein.lstrip('{').rstrip('}')
            treatment = treatmentGroupDict.get(sample)
            labStat = statDict.get(treatment)
            if stripProtein.startswith('MGYG') and labStat != 'Unlabeled':
                splitProtein = stripProtein.split(',')[0]
                countData.append([splitProtein, treatment])
                if enrich <= 100:
                    enrichData.append([splitProtein, treatment, enrich])

        enrichDataDf = pd.DataFrame(enrichData).rename(columns = {0: 'Protein', 1: 'Sample', 2: 'Enrichment'})
        ### Get average enrichment of each protein
        t_avgEnrichmentDf = pd.pivot_table(enrichDataDf, index = 'Protein', columns = 'Sample', values = ['Enrichment'], aggfunc = 'mean').fillna(0)

        abundDataDf = pd.DataFrame(countData).rename(columns = {0: 'Protein', 1: 'Sample'})
        ### Assume each row in SIPdf represents 1 PSM
        spectralCountDf = pd.DataFrame(abundDataDf.groupby(['Sample', 'Protein']).size())
        ### Get sum of spectral counts for each protein
        t_spectralCountDf = pd.pivot_table(spectralCountDf, index = 'Protein', columns = 'Sample', aggfunc = 'sum').fillna(0)
        return t_spectralCountDf, t_avgEnrichmentDf

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
            lineageDict[isolate] = splitLineage[5]
        return lineageDict

    def insertTaxonomy(self, df, dataType, taxonomyDict):
        """
        Insert taxon names into average enrichment and spectral count DataFrames

        Parameters
        ----------
        df : pandas.DataFrame 
            Average enrichment or spectral count dataframe where columns are samples and indices are protein IDs

        dataType: str
            Indicates if aggfunc used to summarize data at the taxon level should be 
            mean (average enrichement) or sum (spectral counts)

        taxonomyDict : dict
            Lookup dict with taxon IDs --> genus names ouput by parseMGYGData()

        Returns 
        -------
        gb : pandas.DataFrame
            DataFrame with average enrichment or spectral count values summarized at the taxon level
        """
        taxonomyData = {}
        for proteinID, sampleData in df.iterrows():
            taxonID = proteinID.split('_')[0]
            taxon = taxonomyDict.get(taxonID)
            ### Need to generate a dictionary that contains protein IDs --> genus names
            taxonomyData[proteinID] = taxon 
        df = df.reset_index()
        df['Taxon'] = df['Protein'].map(taxonomyData)
        df = df.set_index('Taxon').drop('Protein', axis = 1, level = 0).reset_index()
        if dataType == 'spectral_counts':
            gb = df.groupby('Taxon').sum()
        if dataType == 'average_enrichment':
            gb = df.groupby('Taxon').mean()
        return gb

class plotGenera():
    """
    Plot parsed SIP data summarized at the taxon level as a bubble plot

    Attributes:
        averageEnrichmentData (pandas.DataFrame):
            DataFrame with average enrichment values where index is genus name and columns represent samples

        spectralCountData (pandas.DataFrame):
            DataFrame with sum spectral count values where index is genus name and columns represent samples

        sigLabeledTaxa (pandas.DataFrame):
            DataFrame with one column containing names of genera deemed significant based on null distributions

        colormap (LinearSegmentedColormap object):
            Custom colormap for visualization
    """

    def __init__(self, averageEnrichmentData, spectralCountData, sigLabeledTaxa, colormap):
        self.spectralCountData = spectralCountData
        self.averageEnrichmentData = averageEnrichmentData
        self.sigLabeledTaxa = sigLabeledTaxa
        self.colormap = colormap

    def chooseGenera(self):
        """
        Of all significantly labeled taxa, determine which ones should be included in visualization based on thresholds
        
        Parameters
        ----------
        self.averageEnrichmentData : pandas.DataFrame
            DataFrame with average enrichment values where index is genus name and columns represent samples

        self.spectralCountData : pandas.DataFrame
            DataFrame with sum spectral count values where index is genus name and columns represent samples

        Returns
        -------
        taxa2use : list
            List of significantly labeled genera whose average enrichment or spectral count values surpass thresholds

        aeValuesList : list
            List of all average enrichment values that will be included in plot.
            Used for determining max and min average enrichment values. 

        scValuesList : list
            List of all spectral count values that will be included in plot. 
            Used for determining max spectral count value to display in legend. 

        Notes
        -----
        [insert threshold justification]
        """
        taxa2use = []
        aeValuesList = []
        scValuesList = []
        ### Set thresholds
        minSpectralCountThresh = 10
        avgEnThresh = 25
        maxAvgEn = 95
        spectCountThresh = 50
        for (taxon_e, *samples_e), (taxon_s, *samples_s) in zip(self.averageEnrichmentData.itertuples(), self.spectralCountData.itertuples()):
            if taxon_e == taxon_s and taxon_e in self.sigLabeledTaxa['Taxon'].values:
                if (not all(specc < minSpectralCountThresh for specc in samples_s) and any(specc > spectCountThresh for specc in samples_s)) or (not all(specc < minSpectralCountThresh for specc in samples_s) and any(maxAvgEn > avgen > avgEnThresh for avgen in samples_e)):
                    taxa2use.append(taxon_e)
                    for num in samples_e:
                        aeValuesList.append(num)
                    for num in samples_s:
                        scValuesList.append(num)
        return taxa2use, aeValuesList, scValuesList

    def plotBubblePlot(self, enrichmentVals, spectCountVals, taxa):
        """
        Plot visualization

        Parameters
        ----------
        enrichmentVals : list
            List of all average enrichment values that will be included in plot.
        
        spectCountVals : list
            List of all spectral count values that will be included in plot. 

        """
        fig, ax = plt.subplots(figsize = (11.5, 8))
        plt.subplots_adjust(left=0.26)

        vmin = min(enrichmentVals)
        vmax = max(enrichmentVals)
        norm = Normalize(vmin=vmin, vmax=vmax) 
        
        plotTaxa = pd.Series(taxa).drop_duplicates().values.tolist()
        self.spectralCountData.columns = self.spectralCountData.columns.droplevel(0)
        self.averageEnrichmentData.columns = self.averageEnrichmentData.columns.droplevel(0)
        ### Organize columns based on time point so x axis is chronological
        self.spectralCountData = self.spectralCountData.reindex(['C6', 'C12', 'C18', 'C24'], axis=1)
        self.averageEnrichmentData = self.averageEnrichmentData.reindex(['C6', 'C12', 'C18', 'C24'], axis=1)

        ylocs = []
        ylabs = []
        cmin = 0
        cmax = max(spectCountVals)
        minS, maxS = 30, 900 
        for i1, (taxon) in enumerate(plotTaxa):
            # get series with average enrichment for each taxon
            enrichmentData = self.averageEnrichmentData.loc[taxon] 
            # get series with sum spectral counts for each taxon
            countData = self.spectralCountData.loc[taxon] 
            ylocs.append(i1) 
            ylabs.append(taxon)
            ### iterate over each series with enrichment or spectral count data of all samples
            for i2, (eData, cData) in enumerate(zip(enrichmentData.values, countData.values)): 
                ### scale size factor
                sval = minS + (maxS - minS) * (cData - cmin) / (cmax - cmin) if cmax > cmin else minS
                if cData == 0:
                    plt.scatter(i2, i1, s = sval,  c = 'white', norm = norm)
                else:
                    plt.scatter(i2, i1, s = sval, c = [eData], cmap = self.colormap, norm = norm)

        xLocs = []
        xLabs = []
        for idx, sname in enumerate(list(self.spectralCountData.columns)):
            xLocs.append(idx)
            xLabs.append(sname)
        plt.yticks(ylocs, ylabs, fontsize = 16)
        plt.xticks(xLocs, xLabs, fontsize = 16)

        ### Get range of values for spectral count legend
        exampleCounts = [1,round(max(spectCountVals)/2, 0),max(spectCountVals)]
        handles = []
        for c in exampleCounts:
            sval = minS + (maxS - minS) * (c - cmin) / (cmax - cmin) if cmax > cmin else minS
            handles.append(ax.scatter([], [], s=sval, color="gray", alpha=0.6))
        labels = [f"{int(c)}" for c in exampleCounts]
        sm = ScalarMappable(norm=norm, cmap=self.colormap)
        sm.set_array([])
        ax.legend(handles, labels, title="Spectral Count", scatterpoints=1, frameon=False, labelspacing=2.5, bbox_to_anchor=(1.4,1.02))
        fig.colorbar(sm, ax=ax, label='Average Enrichment')
        plt.show()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inFile')
    parser.add_argument('-n', '--namesDict')
    parser.add_argument('-m', '--metadata')
    parser.add_argument('-t', '--sigtaxa')
    args = parser.parse_args()

    colormap = LinearSegmentedColormap.from_list("colorList", ["gold", "orange", "crimson"])
    sipData = pd.read_csv(args.inFile, sep = '\t', header = 0, usecols = [0, 17, 18, 25, 26])
    sigLabeledTaxa = pd.read_csv(args.sigtaxa, names=['Taxon'])

    dataParser = parseSIPData(sipData)
    groupDict, statusDict = dataParser.sampleMetadata(args.namesDict)
    taxonomyLookupDict = dataParser.parseMGYGData(args.metadata)
    spectralCounts, averageEnrichment = dataParser.parseInFile(groupDict, statusDict)

    taxonSpectralCounts = dataParser.insertTaxonomy(spectralCounts, 'spectral_counts', taxonomyLookupDict)
    taxonAverageEnrichment = dataParser.insertTaxonomy(averageEnrichment, 'average_enrichment', taxonomyLookupDict)

    chooser = plotGenera(taxonAverageEnrichment, taxonSpectralCounts, sigLabeledTaxa, colormap)
    genera, avgEnValues, spectCountValues = chooser.chooseGenera()
    chooser.plotBubblePlot(avgEnValues, spectCountValues, genera)


if __name__ == "__main__":
    main()
