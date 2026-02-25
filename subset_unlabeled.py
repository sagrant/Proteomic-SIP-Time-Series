#!/usr/bin/env python3
"""
subset_unlabeled.py 

Purpose: 
    Subset sample-specific Sipros outputs to only include unlabeled PSMs.
    Combine into one DataFrame containing all unlabeled PSMs for all samples.

Inputs:
    - Path to sample-specific Sipros output files

Outputs:
    - DataFrame containing all unlabeled PSMs for all samples.

Usage:
    python subset_unlabeled.py \
        -p [path] \
        -o [output file name]

Notes:
    - Any PSM enriched <2 % is unlabeled
    - Resulting DataFrame will have the same strucutre as Percolator output
"""

import pandas as pd
import numpy as np
import argparse
import os

class subsetUnlabeledProteins():
        """
        Subset concatenated DataFrame so it only contains unlabeled PSMs.

        Attributes:
            sipDfsList (list) : list of DataFrames with unlabeled and labeled PSMs
            There is one DataFrame per sample
        """
        
        def __init__(self, sipDfsList):
            self.sipDfsList = sipDfsList        

        def unlabeledProteinIDs(self):
            """
            Subset dataframes

            Returns
            -------
            unlabeledConcatDf : pandas.DataFrame
                DataFrame containing all unlabeled PSMs for all samples
            """
            concatDf = pd.concat(self.sipDfsList)            
            unlabeledPSMsList = []
            for fields in concatDf.itertuples(index = False):
                if fields[17] >= 2 and fields[18] >= 2:
                    unlabeledPSMsList.append(fields[0])
            concatDf = concatDf.set_index('PSMId')
            unlabeledConcatDf = concatDf.loc[unlabeledPSMsList]
            return unlabeledConcatDf

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path')
    parser.add_argument('-o', '--outFile')
    args = parser.parse_args()

    path = args.path
    pathList = os.listdir(path)
    dfsList = []
    for fname in pathList:
        df = pd.read_csv(f'{path}/{fname}', sep = '\t')
        sampleFile = fname.split('_filtered_psms.tsv')[0]
        df['FileName'] = sampleFile
        dfsList.append(df)

    subsetUnlabeled = subsetUnlabeledProteins(dfsList)
    unlabeledDf = subsetUnlabeled.unlabeledProteinIDs()
    unlabeledDf.to_csv(args.outFile, sep = '\t')

if __name__ == "__main__":
    main()