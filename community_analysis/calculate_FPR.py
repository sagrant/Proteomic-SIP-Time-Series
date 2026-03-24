#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse

'''
calculate_FPR.py

Purpose:
    Calculate False Positive Rate (FPR) for dataset

Inputs:
    - Percolator output (TSV)
    - File names --> sample names lookup table (CSV)
    - Table with number of spectra aquired for each sample, generated with the Aerith package in R (CSV)

Outputs:
    - Dataset FPR printed to STDOUT

Usage:
    python calculate_FDR.py \
        -n [sample names lookup table] \
        -p [percolator output] \
        -s [spectra table]
'''

class calculateDataSetFPR():
    """
    Calculate FPR

    This class computes the dataset FPR based on the number of labeled PSMs
    detected in unlabeled samples that passed Percolator filtering and the total spectra
    acquired for unlabeled samples
    """

    def __init__(self, SIPPsms, sampleSpectra):
        self.SIPPsms = SIPPsms
        self.sampleSpectra = sampleSpectra

    def sampleMetadata(self, namesDictIn):
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
        sampleLookup = pd.read_csv(namesDictIn)
        sampleLookupDict = sampleLookup.to_dict(orient = 'index')
        sampDict = {}
        sampStatDict = {}
        for sampleID, sampleName in sampleLookupDict.items():
            sampDict[sampleName['FileName']] = sampleName['SampleName']
            if '0' in sampleName['SampleName'].split('.')[0]:
                sampStatDict[sampleName['SampleName']] = 'Unlabeled'
            else:
                sampStatDict[sampleName['SampleName']] = 'Labeled'
        return sampDict, sampStatDict

    def PSMCounts(self, sDict, ssDict):
        """
        Count number of labeled PSMs in unlabeled samples 

        Parameters 
        ----------
        sDict : dict
            Dictionary to convert MS/MS file name to sample name
        ssDict : dict
            Dictionary to convert sample name to label status (Labeled/Unlabeled)

        Returns
        -------
        totLabeled_unlabeledSamples : list
            List with counts of labeled PSMs in all unlabeled samples
        """
        gbSIPdf = self.SIPPsms.groupby(self.SIPPsms.iloc[:,2])
        totLabeled_unlabeledSamples = []
        for x, sampleGroup in gbSIPdf:
            sampleName = sDict.get(x)
            status = ssDict.get(sampleName)
            if status == 'Unlabeled':
                totLabeled_unlabeledSamples.append(len(sampleGroup))
        return totLabeled_unlabeledSamples

    def SpectraCounts(self, ssDict):
        """
        Count number of acquired spectra in unlabeled samples 

        Parameters
        ----------
        ssDict : Dict
            Dictionary to convert MS/MS file name to sample name

        Returns
        -------
        totSpectra_unlabeledSamples : int
            Sum of all spectra acquired across all unlabeled samples
        """
        self.sampleSpectra['Status'] = self.sampleSpectra.iloc[:,0].map(ssDict)
        unlabeledSampleSpectra = self.sampleSpectra.set_index('Status').loc['Unlabeled']
        totSpectra_unlabeledSamples = unlabeledSampleSpectra.iloc[:,1].sum()
        return totSpectra_unlabeledSamples

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--namesDict')
    parser.add_argument('-p', '--percolatorOutput')
    parser.add_argument('-s', '--spectra')
    args = parser.parse_args()

    SIPdf = pd.read_csv(args.percolatorOutput, sep = '\t', header = 0, usecols = [0, 25, 26])
    spectraDf = pd.read_csv(args.spectra).dropna()

    calcFPR = calculateDataSetFPR(SIPdf, spectraDf)
    sampleDict, statusDict = calcFPR.sampleMetadata(args.namesDict)
    totalLabeled_unlabeledSamples = calcFPR.PSMCounts(sampleDict, statusDict)
    totalSpectra_unlabeledSamples = calcFPR.SpectraCounts(statusDict)

    FPR = sum(totalLabeled_unlabeledSamples) / totalSpectra_unlabeledSamples
    print('------------------------------------------')
    print(f'FPR for this dataset = {FPR}')
    print('------------------------------------------')
    print('------------------------------------------')

if __name__ == "__main__":
    main()