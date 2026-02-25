#!/usr/bin/env python3
"""
proportion_labeled_PSMs.py

Purpose:
    Visualize and record the proportion of labeled PSMs out of all detected PSMs.

Inputs:
    - Path to directory with sample-specific Sipros output files in TSV format
    - File names --> sample names lookup table (CSV)

Outputs:
    - Bar chart depicting changes in proportion of labeled PSMs over time
    - CSV containing proportion of labeled and unlabeled PSMs, as well as total counts of PSMs

Usage:
    python proportion_labeled_PSMs.py \
        -p [file path] \
        -n [sample names lookup table] \
        -o [output file]

"""

import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import os

class calculateProportions():
    """
    Calculate the proportion of unlabeled and labeled PSMs.
    Record total count of PSMs. 

    Attributes:
        dfsList (list) : List of pandas DataFrames, one per sample
    """

    def __init__(self, dfsList):
        self.dfsList = dfsList

    def sampleMetadata(self, namesDictIn):
        """
        Generate lookup dictionaries containing metadata for all samples

        Parameters
        ----------
        namesDictIn : file path to CSV
            CSV file contains lookup table to convert MS/MS file name --> experimental sample name 

        Returns 
        -------
        groupDict : dict
            Dictionary to convert MS/MS file name to treatment group ID

        orderDict : dict
            Dictionary used to map treatment group order to DataFrame
            Allows for chronological ordering of treatement group-specific values
        
        Notes
        -----
        Samples from this experiment are named like [Sample type][Time point].[Replicate]    
        """
        sampleLookup = pd.read_csv(namesDictIn)
        sampleLookupDict = sampleLookup.to_dict(orient = 'index')
        groupDict = {}
        orderDict = {}
        for sampleID, sampleName in sampleLookupDict.items():
            group = sampleName['SampleName'].split('.')
            groupDict[sampleName['FileName']] = group[0]
            ### Generate some dummy values to associate treatment group names with their chronological order
            if 'C' in sampleName['SampleName']:
                orderDict[group[0]] = int(group[0][1::]) + 1
            if 'S' in sampleName['SampleName']:
                orderDict[group[0]] = (int(group[0][1::])+ 25)*2
        return groupDict, orderDict

    def parseSamples(self, treatmentGroupDict, ordDict):
        """
        Parse sample-specific dataframes and record counts, then calculate proportions

        Parameters
        ----------
        treatmentGroupDict : dict
            Dictionary to convert MS/MS file name to treatment group ID

        ordDict : dict 
            Dictionary used to map treatment group order to DataFrame

        Returns
        -------
        countDf : pandas.DataFrame
            Contains proportions of labeled and unlabeled PSMs, as well as total count
            Summarized at the treatment group level 
        """
        countDataBySample = []
        for sampleData in self.dfsList:
            totalPSMs = 0
            labeledPSMs = 0
            unlabeledPSMs = 0
            for psm, ms1, ms2, protein in sampleData.itertuples(index = False):
                sID = psm.split('.')[0]
                sampleName = treatmentGroupDict.get(sID)
                totalPSMs += 1
                if protein.startswith('{MGYG'):
                    if ms2 >= 2 and ms1 >= 2:
                        labeledPSMs += 1
                    else:
                        unlabeledPSMs += 1
            countDataBySample.append([sampleName, totalPSMs, labeledPSMs, unlabeledPSMs])

        countDf = pd.DataFrame(countDataBySample, columns = ['Group', 'Total_PSMs', 'Labeled_PSMs', 'Unlabeled_PSMs'])
        countDf['Order'] = countDf['Group'].map(ordDict)
        countDf = countDf.sort_values(by = 'Order')
        countDf = countDf.drop(['Order'], axis = 1)
        countDf['Proportion_Labeled'] = (countDf['Labeled_PSMs'] / countDf['Total_PSMs'])
        countDf['Proportion_Unlabeled'] = (countDf['Unlabeled_PSMs'] / countDf['Total_PSMs'])
        return countDf

    def plotProportions(self, sampleData):
        """
        Plot proportion of labeled PSMs out of total detected PSMs

        Parameters:
        ----------
        sampleData : pandas.DataFrame
            Contains proportions of labeled and unlabeled PSMs, as well as total count
        """
        gbGroup = sampleData.groupby('Group', sort = False).mean().round(decimals=3)

        fig, axs = plt.subplots()

        labSTDVs = []
        for groupData in sampleData.groupby('Group', sort = False):
            stdv_lab = np.std(groupData[1]['Labeled_PSMs']/groupData[1]['Total_PSMs'])
            labSTDVs.append(stdv_lab)

        axs.bar(list(range(len(gbGroup.index))), gbGroup['Proportion_Labeled'], yerr = labSTDVs, color = 'slategrey')
        axs.set_xticks(list(range(len(gbGroup.index))))
        axs.set_xticklabels(gbGroup.index)
        axs.set_title("Mean Proportion of Labeled PSMs")
        axs.set_xlabel("Treatment Group")
        axs.set_ylabel("Mean Proportion")

        plt.tight_layout()
        plt.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path')
    parser.add_argument('-n', '--names')
    parser.add_argument('-o', '--outFile')
    args = parser.parse_args()

    path = args.path
    pathList = os.listdir(path)
    dfsList = []
    for fname in pathList:
        df = pd.read_csv(f'{path}/{fname}', sep = '\t', usecols = [0, 17, 18, 26])
        dfsList.append(df)

    calc = calculateProportions(dfsList)
    gDict, oDict = calc.sampleMetadata(args.names)
    counts = calc.parseSamples(gDict, oDict)
    calc.plotProportions(counts)
    counts.to_csv(args.outFile)

if __name__ == "__main__":
    main()