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
import random
import itertools
import matplotlib.pyplot as plt

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
    orderDict = {}
    for sampleID, sampleName in sampleLookupDict.items():
        sampleNameDict[sampleName['FileName']] = sampleName['SampleName']
        orderDict[sampleName['SampleName']] = int(sampleName['SampleName'].split('.')[1])
    return sampleNameDict, orderDict

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

    def parse13CEnrichment(self, n, oDict):
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
        random.seed(42)
        plotEnrichmentData = []
        for enrichmentData in self.dfsList:
            graphName = []
            enrichmentValues = []
            subsampleEnrichData = enrichmentData.sample(n=n, replace=False, random_state=42)
            for psm, enrichment1, enrichment2, prot, sample in subsampleEnrichData.itertuples(index = False):
                group = sample.split('.')[0]
                mouse = int(sample.split('.')[1])
                ### Only save enrichment values of mouse proteins
                if prot.startswith('{MGYG'):
                    if enrichment2 >= 2 and enrichment1 >= 2 and enrichment1 <= 100:
                        enrichmentValues.append(enrichment1)
                        graphName.append(sample)
            plotEnrichmentData.append([enrichmentValues, graphName[0], group, mouse])
        df = pd.DataFrame(plotEnrichmentData, columns = ['Values', 'Sample', 'Group', 'Mouse'])
        df['Order'] = df['Sample'].map(oDict)
        ### Sort samples in chronological order
        df = df.sort_values(by = 'Order')
        return df

    def generateColorMap(self, mouseList):
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
            "yellowgreen",
            "forestgreen",
            "blue",
            "firebrick",
            "darkgoldenrod",
            "darkorange",
            "orangered",
            "darkmagenta",
            "maroon",
            "saddlebrown",
            "palevioletred",
            "teal",
            "magenta",
            "thistle"]
        mouseColMap = {}
        for mouse, color in zip(mouseList, colors):
            mouseColMap[mouse] = color
        return mouseColMap
    
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
        plt.subplots_adjust(left = 0.09, bottom = 0.06, top = 0.95, right = 0.85)

        enrichData['Color'] = enrichData['Mouse'].map(colorDict)
        ### plot each group, each with 3 mice, on its own subplot
        gbGroup = enrichData.groupby('Group', sort = False)

        for (groupName, enrichmentValueData), axis in zip(gbGroup, ax.flatten()):
            for values, samp, gp, mouse, order, colr in enrichmentValueData.itertuples(index = False):
                axis.hist(values, bins = 50, alpha = 0.9, log = True, label = mouse, color = colr, histtype = 'step', linewidth = 1.3)
                axis.set_title(gp)
                axis.set_xticks([0, 20, 40, 60, 80, 100])
                axis.set_xlim(0,103)
                axis.set_ylim(1, 150)
        
        ax[4].set_xlabel('Enrichment (%)')
        ax[2].set_ylabel(r'$\mathrm{Log}_{10}$ PSM Count')

        ### Generate legend so there are no duplicate entries 
        ### and mice are listed in chronological order
        legendInfo = []
        legendOrder = {}
        for i, a in enumerate(ax.ravel()):
            h, l = a.get_legend_handles_labels()
            l.sort()
            for hand, lab in zip(h, l):
                if lab not in legendInfo:
                    legendOrder[lab] = i
                    legendInfo.append([hand, lab])
        legendDf = pd.DataFrame(legendInfo, columns = ['Handles', 'Labels']).drop_duplicates('Labels')
        legendDf['Order'] = legendDf['Labels'].map(legendOrder)
        legendDf = legendDf.sort_values(by = 'Order')
        handles = legendDf['Handles'].values.tolist()
        labels = legendDf['Labels'].values.tolist()

        plt.legend(handles, labels, bbox_to_anchor=(1.19,4), title='Replicate', fontsize = 10)
        plt.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path')
    parser.add_argument('-n', '--names')
    args = parser.parse_args()

    sDict, orDict = sampleMetadata(args.names)

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
    data = community13C.parse13CEnrichment(subsampleValue, orDict)
    colorProfile = community13C.generateColorMap(data['Mouse'].unique())
    community13C.plotEnrichmentDistributions(data, colorProfile)
    
if __name__ == "__main__":
    main()