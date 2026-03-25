#!/usr/bin/env python3
"""

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
    sDict = {}
    sampStatDict = {}
    for sampleID, sampleName in sampleLookupDict.items():
        group = sampleName['SampleName'].split('.')
        sDict[sampleName['FileName']] = sampleName['SampleName']
        if '0' in sampleName['SampleName'].split('.')[0]:
            sampStatDict[group[0]] = 'Unlabeled'
        else:
            sampStatDict[group[0]] = 'Labeled'
    return sDict, sampStatDict


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

    def __init__(self, siprosData):
        self.siprosData = siprosData

    def computeSpectralCounts(dataList, colName):
        abundanceData = pd.DataFrame(dataList).rename(columns = {0:'Genus'})
        spectralCountDf = abundanceData.groupby('Genus').size().reset_index().rename(columns = {0:colName})
        return spectralCountDf

    def parseSiprosData(self, taxonomyDict, k0Dict, keggFunctDict):
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
                if k0function:
                    funct = keggFunctDict.get(k0function['kegg'], 'Absent')
                    k0taxonData.append([taxon, funct])
                    sampleData.append([funct.replace(' ', '\n'), sample])
                if not k0function:
                    k0taxonData.append([taxon, 'Absent'])
                if ms2 >= 2:
                    enrichmentData.append([taxon, ms2])
                    labeledSpectralCountData.append(taxon)
        
        k0taxonDf = pd.DataFrame(k0taxonData, columns=['Genus', 'KEGG_Function']).fillna(0)
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

    sampDict, statDict = sampleMetadata(args.namesDict)
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
        sampleStat = statDict.get(sid)
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
