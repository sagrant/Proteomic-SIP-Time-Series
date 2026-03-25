#!/usr/bin/env python3
"""
trait_data.py

Purpose:
    Generate dataframe required to make phylogenetic tree heatmap

Inputs:
    - Path to directory with Sipros output files
    - File names --> sample names lookup table (CSV)
    - MGYG proteome metadata with taxon IDs --> taxon lineage (CSV)
    - DataFrame containing all MGYG protein IDs associated with enzymes of interest (CSV)

Outputs:
    - Trait data to be used for make phylogenetic tree heatmap

Usage:
    python generate_phyloseq.py \
        -p [file path] \
        -n [sample dictionary] \
        -m [MGYG metadata] \
        -k [KEGG functions metadata]
        -o [output file]
"""
import pandas as pd
import numpy as np
import argparse
import matplotlib.pyplot as plt
import os
import seaborn as sns


def sampleMetadata(namesDictIn):
    """
    Generate lookup dictionaries containing metadata for all samples

    Parameters
    ----------
    namesDictIn : file path
        CSV file contains lookup table to convert MS/MS file name --> experimental sample name 

    Returns 
    -------
    sDict : dict
        Dictionary to convert MS/MS file name to sample name

    Notes
    -----
    Samples from this experiment are named like [Sample type][Time point].[Replicate]    
    """
    sampleLookup = pd.read_csv(namesDictIn)
    sampleLookupDict = sampleLookup.to_dict(orient = 'index')
    sDict = {}
    for sampleID, sampleName in sampleLookupDict.items():
        sDict[sampleName['FileName']] = sampleName['SampleName']
    return sDict


def parseMGYGData(metadata):
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


class generateTraitData():
    """
    Generate trait data 

    Resulting dataframe contains labeled spectral counts, total spectral counts,
    average enrichment, and spectral counts of enzymes of interest for all detected genera

    Attributes:
        siprosData (pandas.DataFrame) : Concatenated sipros output for all samples
    """

    def __init__(self, siprosData):
        self.siprosData = siprosData

    def computeSpectralCounts(dataList, colName):
        """
        Calculate total spectral counts or labeled spectral counts of all detected genera

        Parameters
        ----------
        dataList : list
            List of genera names, each item in list represents 1 spectral count

        colName : str
            Name of spectral count column in resulting dataframe 

        Returns
        -------
        spectralCountDf : pandas.DataFrame
            Dataframe with either total spectral counts or labeled spectral counts of all 
            detected genera
        """
        abundanceData = pd.DataFrame(dataList).rename(columns = {0:'Genus'})
        spectralCountDf = abundanceData.groupby('Genus').size().reset_index().rename(columns = {0:colName})
        return spectralCountDf

    def parseSiprosData(self, taxonomyDict, k0Dict, keggFunctDict):
        """
        Parse total proteome and save functional information, spectral counts, and average
        enrichment values

        Parameters
        ----------
        taxonomyDict : dictionary
            Dictionary used to retrieve genus names, where keys are MGYG taxon ID and 
            values are genus names

        k0Dict : dictionary 
            Dictionary used to retreive Kegg IDs for enzymes of interest. Keys are MGYG ID
            and values are K0 number

        keggFunctDict : dictionary 
            Hard-coded dictionary where keys are K0 number and values are the functional name
            of the enzyme 

        Returns
        -------
        mergedTraitData : pandas.DataFrame
            Dataframe containing all data that will be used to make phylogenetic tree heatmap
        """
        spectralCountData = []
        labeledSpectralCountData = []
        k0taxonData = []
        enrichmentData = []
        sampleData = []
        for psm, ms2, protein, sample in self.siprosData.itertuples(index = False):
            stripProtein = protein.lstrip('{').rstrip('}')
            if stripProtein.startswith('MGYG'):
                splitProtein = stripProtein.split(',')[0]
                taxonID = splitProtein.split('_')[0]
                taxon = taxonomyDict.get(taxonID)
                spectralCountData.append(taxon)
                k0function = k0Dict.get(splitProtein)
                ### If enzyme is associated with any gene in current taxon's genome, save it 
                if k0function:
                    funct = keggFunctDict.get(k0function['kegg'])
                    k0taxonData.append([taxon, funct])
                    ### Save function and sample names to visualize expression of these enzymes over time
                    sampleData.append([funct.replace(' ', '\n'), sample])
                ### If enzyme is NOT associated with genes encoded by current taxon, save placeholder string
                if not k0function:
                    k0taxonData.append([taxon, 'Absent'])
                if ms2 >= 2:
                    enrichmentData.append([taxon, ms2])
                    labeledSpectralCountData.append(taxon)
        
        k0taxonDf = pd.DataFrame(k0taxonData, columns=['Genus', 'KEGG_Function']).fillna(0)
        ### Use size because we assume each row = 1 spectral count
        gbk0Taxon = k0taxonDf.groupby(['Genus', 'KEGG_Function']).size().reset_index().rename(columns = {0: 'Spectral_Count'})
        pivotk0Tax = pd.pivot_table(gbk0Taxon, index = 'KEGG_Function', columns= 'Genus', values = 'Spectral_Count', aggfunc=sum).T.fillna(0).reset_index()
        
        enrichDataDf = pd.DataFrame(enrichmentData).rename(columns = {0: 'Genus', 1: 'Enrichment'})
        avgEnrichmentDf = enrichDataDf.groupby('Genus').mean()
        
        totalSpectralCountDf = generateTraitData.computeSpectralCounts(spectralCountData, 'Total_Spectral_Count')    
        spectralCountDf = generateTraitData.computeSpectralCounts(labeledSpectralCountData, 'Labeled_Spectral_Count')    

        mergeLabeled = spectralCountDf.merge(avgEnrichmentDf.reset_index(), on = 'Genus',  how = 'outer')
        mergeTotals = mergeLabeled.merge(totalSpectralCountDf, on = 'Genus', how = 'outer').fillna(0)

        mergedTraitData = pivotk0Tax.merge(mergeTotals, on = 'Genus', how = 'outer').fillna(0)
        return mergedTraitData
    
            
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--path')
    parser.add_argument('-n', '--namesDict')
    parser.add_argument('-m', '--metadata')
    parser.add_argument('-k', '--keggDict')
    parser.add_argument('-o', '--outFile')
    args = parser.parse_args()

    sampDict = sampleMetadata(args.namesDict)
    lingDict = parseMGYGData(args.metadata)

    k0sDf = pd.read_csv(args.keggDict, sep = ',').set_index('gene_id')
    k0sLookupDict = k0sDf.to_dict(orient = 'index')
    k0sFunctionDict = {'K00929': 'butyrate kinase', 'K00248':'butyryl-CoA dehydrogenase', 'K00074':'3-hydroxybutyryl-CoA dehydrogenase'}

    path = args.path
    pathList = os.listdir(path)
    dfs = []
    for fHandle in pathList:
        sid = fHandle.split('_filtered_psms.tsv')[0]
        sample = sampDict.get(sid)
        ### Only include cecum samples
        if 'C' in sample:
            df = pd.read_csv(f'{path}/{fHandle}', sep = '\t', usecols = [0, 18, 26])
            df['Sample'] = sample
            dfs.append(df)

    concatDf = pd.concat(dfs)

    traitDataGenerate = generateTraitData(concatDf)
    outData = traitDataGenerate.parseSiprosData(lingDict, k0sLookupDict, k0sFunctionDict)
    outData.to_csv(args.outFile)

if __name__ == "__main__":
    main()
