#!/usr/bin/env python3
"""
enrichment_distributions.py

Purpose:
    Plot 13C enrichment distributions of all labeled PSMs to study
    labeling of the community over time

Inputs:
    - Path to directory with sample-specific Sipros output files
    - File names --> sample names lookup table (CSV)

Outputs:
    - Enrichment distribution visualization

Usage:
    python enrichment_distributions.py \
        -p [file path] \
        -n [sample lookup table]

"""
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

def sampleMetadata(namesDict):
    """
    Generate lookup dictionaries containing metadata for all samples

    Parameters
    ----------
    namesDictIn : file path to CSV
        CSV file contains lookup table to convert MS/MS file name --> experimental sample name 

    Returns 
    -------
    sampleNameDict : dict
        Dictionary to convert MS/MS file name to sample name
    
    sampStatDict : dict
        Dictionary that can be used to get label status of each sample (labeled/unlabeled)

    Notes
    -----
    Samples from this experiment are named like [Sample type][Time point].[Replicate]    
    """
    sampleLookup = pd.read_csv(namesDict)
    sampleLookupDict = sampleLookup.to_dict(orient = 'index')
    sampleNameDict = {}
    for sampleID, sampleName in sampleLookupDict.items():
        sampleNameDict[sampleName['FileName']] = sampleName['SampleName']
    return sampleNameDict

class communityEnrichment():
    """
    Parse data and plot enrichment distributions of all labeled PSMs

    This class parses the Sipros outputs and plots 13C enrichment values 
    detected in PSMs from all cecum samples as histograms
    """

    def __init__(self, dfsList):
        self.dfsList = dfsList

    def calculateSubsampleSize(self, lineCounts):
        """
        Subsample number of PSMs in dataframes

        Parameters
        ----------
        lineCounts : list
            Nested list containing sample name and the number of PSMs detected 
            in that sample

        Returns
        -------
        subsampleSize : int
            number that all DataFrames should be subsampled to

        Notes
        -----
        To facilitate a fair comparison across samples, it is necessary to
        account for their variable sampling depths. The number of PSMs incorporated
        into this analysis are subsampled to the number of PSMs detected in the 
        sample with fewest total PSMs
        """
        sampleSizesDf = pd.DataFrame(lineCounts, columns = ['Sample', 'Sample_Size'])
        subsampleSize = sampleSizesDf['Sample_Size'].min()
        return subsampleSize

    def parse13CEnrichment(self, n):
        """
        Parse enrichment values of labeled PSMs for all samples

        Parameters
        ----------
        n : int
            Value to subsample all dataframes to 

        oDict : dict
            Dictionary with keys as sample names and values 
            as integers. Used to organize the samples in 
            chronological order

        Returns
        -------
        df : pandas.DataFrame
            DataFrame containing all enrichment values of labeled PSMs, 
            corresponding sample, corresponding treatment group, and 
            mouse ID
        """
        plotEnrichmentData = []
        for enrichmentData in self.dfsList:
            graphName = []
            enrichmentValues = []
            subsampleEnrichData = enrichmentData.sample(n=n, replace=False, random_state=42)
            for psm, enrichment1, enrichment2, prot, sample in subsampleEnrichData.itertuples(index = False):
                group = sample.split('.')[0]
                mouse = int(sample.split('.')[1])
                ### Only save enrichment values of microbial proteins
                if prot.startswith('{MGYG'):
                    if enrichment2 >= 2 and enrichment1 >= 2 and enrichment1 <= 100:
                        enrichmentValues.append(enrichment2)
                        graphName.append(sample)
            plotEnrichmentData.append([enrichmentValues, graphName[0], group, mouse])
        df = pd.DataFrame(plotEnrichmentData, columns = ['Values', 'Sample', 'Group', 'Mouse'])
        ### Sort samples in chronological order
        df = df.sort_values(by = 'Mouse')
        return df

    def generateColorMap(self, groupList):
        """
        Make custom color map so all mice are assoicated with one color on plot
        
        Parameters
        ----------
        mouseList : list
            list of all mouse IDs in dataset

        Returns
        -------
        mouseColMap : dict
            Dictionary with mouse ID as keys and color as values
        """
        colors = [
            "lightblue",
            "slategrey",
            "dodgerblue",
            "teal",
            "navy",]
        colMap = {}
        for timePointGroup, color in zip(groupList, colors):
            colMap[timePointGroup] = color
        return colMap
    
    def plotEnrichmentDistributions(self, enrichData, colorDict):
        """
        Generate histograms

        Parameters
        ----------
        enrichData : pandas.DataFrame 
            DataFrame containing all enrichment values of labeled PSMs, 
            corresponding sample, corresponding treatment group, and 
            mouse ID
        
        colorDict : dict
            Color dict used to map colors to specific mice
        """
        fig, ax = plt.subplots(5, 1, figsize = (7, 8))
        plt.subplots_adjust(left = 0.11, bottom = 0.06, top = 0.95, right = 0.95, hspace=0.5)

        enrichData['Color'] = enrichData['Group'].map(colorDict)
        ### plot each group, each with 3 mice, on its own subplot

        gbGroup = enrichData.groupby('Group', sort = False)

        for (groupName, enrichmentValueData), axis in zip(gbGroup, ax.flatten()):
            concatGroupValues = np.concatenate(enrichmentValueData['Values'].values)
            for values, samp, gp, mouse, colr in enrichmentValueData.itertuples(index = False):
                panelName = f'{gp[1::]} hrs'
                sns.kdeplot(concatGroupValues, ax = axis, color = colr, linewidth=1.3, bw_adjust=0.27, clip=(0, 100))
                axis.set_title(panelName)
                axis.set_xticks([0, 20, 40, 60, 80, 100])
                axis.set_xlim(0,103)
        
        ax[4].set_xlabel('Enrichment (%)')
        plt.show()

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path')
    parser.add_argument('-n', '--names')
    args = parser.parse_args()

    sDict = sampleMetadata(args.names)

    path = args.path
    pathList = os.listdir(path)
    dfs = []
    cLines = []
    for fHandle in pathList:
        sid = fHandle.split('_filtered_psms.tsv')[0]
        sample = sDict.get(sid)
        ### Only include cecum samples
        if 'C' in sample:
            df = pd.read_csv(f'{path}/{fHandle}', sep = '\t', usecols = [0, 17, 18, 26])
            df['Sample'] = sample
            ### Length of the data frame is used to normalize data based on sampling depth 
            cLines.append([sample, len(df)])
            dfs.append(df)
    
    community13C = communityEnrichment(dfs)
    subsampleValue = community13C.calculateSubsampleSize(cLines)
    data = community13C.parse13CEnrichment(subsampleValue)
    colorProfile = community13C.generateColorMap(data['Group'].unique())
    community13C.plotEnrichmentDistributions(data, colorProfile)
    
if __name__ == "__main__":
    main()