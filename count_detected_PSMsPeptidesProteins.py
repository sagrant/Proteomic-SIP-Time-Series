#!/usr/bin/env python3
"""
count_detected_PSMsPeptidesProteins.py

Purpose:
    Count the number of labeled PSMs, peptides, and proteins detected in all samples.

Inputs:
    - Percolator output (TSV)

Outputs:
    - Table with counts of labeled PSMs, peptides, and proteins in CSV format

Usage:
    python count_detected_PSMsPeptidesProteins.py \
        -i [percolator output] \
        -n [sample names lookup table] \
        -o [output counts table]

"""

import pandas as pd
import numpy as np
import argparse

class countDetected():
    """
    Count detected labeled PSMs, peptides, and proteins in dataset

    This class assumes each row in Percolator output represents one PSM
    The number of unique peptides in Peptide column represents total peptides
    The number of unique proteins in Proteins column represents total proteins

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
    

    def countPSMsPeptidesProteins(self, treatGroupDict, treatmentOrderDict):
        """
        Calculate counts 

        Parameters 
        ----------
        treatGroupDict : dict
            Dictionary to convert MS/MS file name to treatment group ID

        treatmentOrderDict : dict 
            Dictionary used to map treatment group order to DataFrame

        Returns
        -------
        labeledCountsDf : pandas.DataFrame
            Contains counts of labeled PSMs, peptides, and proteins detected in each treatment group
            Organized by sample type and in chronological order

        """
        gbSample = self.SIPdf.groupby('FileName')
        LABsampleValues = []
        for sample, data in gbSample:
            tgroupName = treatGroupDict.get(sample)

            labeledPSMs = len(data)
            labeledPeptides = len(data.drop_duplicates('Peptide'))
            labeledProteins = len(data.drop_duplicates('Proteins'))
            LABsampleValues.append([tgroupName, labeledPSMs, labeledPeptides, labeledProteins])

        labeledCountsDf = pd.DataFrame(LABsampleValues, columns = ['TreatmentGroup', 'Labeled PSMs', 'Labeled Peptides', 'Labeled Proteins'])
        labeledCountsDf['Order'] = labeledCountsDf['TreatmentGroup'].map(treatmentOrderDict)
        ### Organize rows in chronological order and sum counts within treatment groups
        labeledCountsDf = labeledCountsDf.sort_values(by = 'Order').drop('Order', axis = 1).groupby('TreatmentGroup', sort=False).sum()
        return labeledCountsDf

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--inFile')
    parser.add_argument('-n', '--sampleNames')
    parser.add_argument('-o', '--outFile')
    args = parser.parse_args()

    sipData = pd.read_csv(args.inFile, sep = '\t', usecols = ['PSMId', 'MS1IsotopicAbundances', 'MS2IsotopicAbundances', 'Peptide', 'Proteins', 'FileName'])
    
    counter = countDetected(sipData)
    groupLookupDict, orderLookupDict = counter.sampleMetadata(args.sampleNames)
    outData = counter.countPSMsPeptidesProteins(groupLookupDict, orderLookupDict)
    outData.to_csv(args.outFile)

if __name__ == "__main__":
    main()