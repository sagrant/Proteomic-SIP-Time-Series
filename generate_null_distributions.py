#!/usr/bin/env python3
"""
generate_null_distributions.py

Purpose:
    Plot null distributions of labeled PSMs for most abundant genera.
    Determine which genera are significantly labeled i.e., their empirical counts of labeled 
    PSMs surpasses null distribution maximum.

Inputs:
    - Path to sample-specific Sipros output files
    - File names --> sample names lookup table (CSV)
    - MGYG proteome metadata with taxon IDs --> taxon lineage 

Outputs:
    - CSV containing names of significantly labeled taxa

Usage:
    python generate_null_distributions.py \
        -p [file path] \
        -n [sample lookup table] \
        -t [MGYG metadata table] \
        -o [output file]

Notes:
    Microbial PSMs are aggregated across all labeled cecum samples and each PSM gets assigned a 
    binary label status (labeled or unlabeled) based on isotopic enrichment threshold of 2%. 
    The total number of observed PSMs was calculated (N) for each genus. To estimate the expected 
    number of labeled PSMs that would be assigned to a genus under random conditions, N PSMs get 
    randomly sampled without replacement from the aggregated dataset, and the number of labeled 
    PSMs in the random sample was recorded. This test gets conducted 1,000 times per genus. Observed 
    labeled PSM counts were compared to the median number of labeled PSMs in the null distribution. 

"""
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import argparse

class nullDistributions():
    """
    Generate and plot null distributions

    This class parses the Sipros outputs, calculates null distributions, and determines
    which genera were significantly labeled
    """

    def sampleMetadata(self, namesDict):
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
        sampStatDict = {}
        for sampleID, sampleName in sampleLookupDict.items():
            sampleNameDict[sampleName['FileName']] = sampleName['SampleName']
            if '0' in sampleName['SampleName'].split('.')[0]:
                sampStatDict[sampleName['FileName']] = 'Unlabeled'
            else:
                sampStatDict[sampleName['FileName']] = 'Labeled'
        return sampleNameDict, sampStatDict


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
            lineageDict[isolate] = splitLineage[5]
        return lineageDict

    def parseSIPfile(self, dfs, linDict):
        """
        Parse list of sample-specific DataFrames containing Sipros outputs

        Parameters
        ----------
        dfs : list 
            list of pandas.DataFrames containing Sipros outputs for all samples

        linDict : dict
            Dictionary to get genus names based on taxon ID

        Returns
        -------
        psmsDf : pandas.DataFrame
            DataFrame containing taxon ID, genus name, and label status
            Label status is a binary variable where 1 = labeled (enriched over 2%) and 0 = unlabeled 
        """
        psmData = []
        for sampleData in dfs:
            for psm, ms1, ms2, protein, sample in sampleData.itertuples(index = False):
                ### Only include microbial proteins
                if protein.startswith('{MGYG'):
                    taxonID = protein.lstrip('{').rstrip('}').split(',')[0].split('_')[0]
                    genus = linDict.get(taxonID)
                    if ms2 >= 2 and ms1 >= 2 and ms1 <= 100:
                        psmData.append([taxonID, genus, 1])
                    else:
                        psmData.append([taxonID, genus, 0])
        psmsDf = pd.DataFrame(psmData, columns=['Taxon_ID', 'Genus', 'Label_Status'])
        return psmsDf

    def chooseTaxa(self, df):
        """
        Identify taxa to plot null distributions for 

        Parameters
        ----------
        df : pandas.DataFrame
            DataFrame returned by parseSIPfile() containing all genera 
            in dataset, taxon ID, and label status 

        Returns
        -------
        abundantGenera : list
            List of most 10 abundant genera in dataset
        """
        gb = df.groupby('Genus').size()
        abundantGenera = gb.nlargest(10).index.tolist()
        return abundantGenera
        
    def calculateNullDist(self, SIPDf, genera):
        """
        Calculate Null Distributions for select genera

        Parameters
        ----------
        SIPDf : pandas.DataFrame
            DataFrame returned by parseSIPfile() containing all genera 
            in dataset, taxon ID, and label status

        genera : list 
            List of genera to calculate distributions for 

        Returns 
        -------
        dists : list
            List of pandas.Series where each series' values is a null distribution 
            and series name is genus

        queryTrue : dict
            Dictionary used to get empirical count of labeled PSMs for each genus
        """
        permutations = 1000
        rng = np.random.default_rng(42)
        labelStatArray = SIPDf['Label_Status'].to_numpy(dtype=np.int8)
        ### n = total number of PSMs in dataset 
        n = labelStatArray.size
        dists = []
        queryTrue = {}
        for taxon in genera:
            ### Take subset of PSMs assigned to current genus 
            taxonData = SIPDf[SIPDf['Genus'] == taxon]
            ### Record empirical count of labeld PSMs
            queryTrue[taxon] = taxonData['Label_Status'].sum()
            ### Initialize array with "garbage" values that will be overwritten
            nullCounts = np.empty(permutations, dtype=np.int16)
            for i in range(permutations):
                ### Take X random indices from 0 - n where X = total number of PSMs assigned to taxon
                ### This determines which PSMs get randomly sampled to generate null distribution
                idx = rng.choice(n, size=len(taxonData), replace=False)
                ### Populate array with sum of labeled PSMs in each randomly sampled null distribution
                nullCounts[i] = labelStatArray[idx].sum()
            nullDistribution = pd.Series(nullCounts, name = taxon)
            dists.append(nullDistribution)
        return dists, queryTrue
    
    def plotDistributions(self, distributions, empirical):
        """
        Plot null distribution histograms for taxa that were significantly labeled

        Parameters
        ----------
        distributions : list
            List of pandas.Series containing null distributions for genera

        empirical : dict
            Dictionary with genus name --> empirical count of labeled PSMs

        Returns
        -------
        sigLabeledSer : pandas.Series
            Series containing names of genera that are significantly labeled
        """
        sigLabeled = []
        for d in distributions:
            empiricalCount = empirical.get(d.name)
            if empiricalCount > max(d):
                plt.hist(d, bins=15, color = 'slategrey')
                plt.title(d.name)
                plt.vlines(empiricalCount, 0, 275, color = 'crimson', label = 'true number labeled')
                plt.vlines(np.median(d), 0, 275, color = 'navy', label = 'true number labeled')
                plt.show()
                sigLabeled.append(d.name)
        sigLabeledSer = pd.Series(sigLabeled)
        return sigLabeledSer

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path')
    parser.add_argument('-n', '--names')
    parser.add_argument('-t', '--MGYGMetadata')
    parser.add_argument('-o', '--outFile')
    args = parser.parse_args()

    nullDists = nullDistributions()
    sDict, sstatDict = nullDists.sampleMetadata(args.names)
    isolatesDict = nullDists.parseMGYGData(args.MGYGMetadata)

    path = args.path
    pathList = os.listdir(path)
    dfs = []
    for fHandle in pathList:
        sid = fHandle.split('_filtered_psms.tsv')[0]
        sample = sDict.get(sid)
        sampleStat = sstatDict.get(sid)
        ### Only include labeled cecum samples
        if sampleStat != 'Unlabeled' and 'C' in sample:
            df = pd.read_csv(f'{path}/{fHandle}', sep = '\t', usecols = [0, 17, 18, 26])
            df['Sample'] = sample
            dfs.append(df)

    sipData = nullDists.parseSIPfile(dfs, isolatesDict)
    taxa = nullDists.chooseTaxa(sipData)
    distList, trueDict = nullDists.calculateNullDist(sipData, taxa)
    sigTaxa = nullDists.plotDistributions(distList, trueDict)
    sigTaxa.to_csv(args.outFile, index = False)

if __name__ == "__main__":
    main()