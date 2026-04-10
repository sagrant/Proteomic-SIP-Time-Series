#!/usr/bin/env python3
import pandas as pd
import os 
import matplotlib.pyplot as plt
import numpy as np
import argparse

"""
collector_curve.py

Purpose:
    Visualize how number of detected genera scales with number of PSMs detected
    Goal is to demonstrate the proteome was adequately sampled to identify all detectable genera

Inputs:
    - Path to directory with sample-specific Sipros output files
    - MGYG proteome metadata with taxon IDs --> taxon lineage 

Usage:
    python collector_curve.py \
        -p [file path] \
        -t [MGYG metadata table] \
"""

class plotCurve():
    """
    Calculate collectors curve and plot results
    """

    def __init__(self, concatDf):
        self.concatDf = concatDf

    def parseMGYGData(self, MGYGMeta):
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

        fullLineageDict : dict
            Lookup dict to get full taxonomic lineage for each taxon ID
        """
        metadataDf = pd.read_csv(MGYGMeta, sep = '\t', header = 0, usecols = ['Genome', 'Lineage'])
        lineageDict = {}
        for isolate, lineage in metadataDf.itertuples(index = False):
            splitLineage = lineage.split(';')
            ### splitLineage is a list that can be indexed to get any taxonomic rank in the lineage
            ### Index 5 corresponds to genera 
            lineageDict[isolate] = splitLineage[5].split('__')[1]
        return lineageDict

    def calcCurve(self, mgygDict):
        """
        Determine how number of detected genera scales with number of PSMs detected
        
        Parameters
        ----------
        mgygDict : dict
            Dictionary to convert protein IDs to genus names

        Returns
        -------
        countDetected : list 
            List that tracks how many taxa were seen with each detected PSM

        countMicrobialPSMs : int
            Integer count of total number of microbial PSMs detected in these data 
        """
        rng = np.random.default_rng(0)
        permuteRows = rng.permutation(self.concatDf.to_numpy())
        seenTaxa = set()
        countDetected = []
        countMicrobialPSMs = 0
        for psmid, prot in permuteRows:
            splitProtein = prot.lstrip('{').rstrip('}').split(',')[0].split('_')[0]
            ### Only count non-degenerate, microbial PSMs 
            if splitProtein.startswith('MGYG'):
                countMicrobialPSMs += 1
                taxon = mgygDict.get(splitProtein)
                ### Save unique taxa that get detected
                seenTaxa.add(taxon)
                ### Record how the number of unique taxa changes with each iteration 
                ### Each iteration represents one PSM
                countDetected.append(len(seenTaxa))
        return countDetected, countMicrobialPSMs
    
    def plotData(self, x, y):
        """
        Plot collectors curve

        Parameters
        ----------
        x : int
            Number of microbial PSMs detected

        y : list
            List of integers that reflects number of unique taxa detected with each 
            detected PSM
        """
        plt.plot(list(range(x)), y, color = 'lightslategrey')
        plt.xlabel('Number of PSMs Detected')
        plt.ylabel('Number of Genera Detected')
        plt.yticks([0, 20, 40, 60, 80], ['0', '20', '40', '60', '80'])
        plt.xticks([0, 200000, 400000, 600000, 800000], ['0', '200000', '400000', '600000', '800000'])
        plt.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path')
    parser.add_argument('-t', '--MGYGMetadata')
    args = parser.parse_args()

    path = args.path
    pathList = os.listdir(path)
    dfsList = []
    for fname in pathList:
        df = pd.read_csv(f'{path}/{fname}', sep = '\t', usecols=[0, 26])
        dfsList.append(df)

    concatDf = pd.concat(dfsList)

    plotter = plotCurve(concatDf)
    lineageDictionary = plotter.parseMGYGData(args.MGYGMetadata)
    taxaCounts, PSMCounts = plotter.calcCurve(lineageDictionary)
    plotter.plotData(PSMCounts, taxaCounts)

if __name__ == "__main__":
    main()