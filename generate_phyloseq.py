
#!/usr/bin/env python3
"""
generate_phyloseq.py

Purpose:
    Step 1 in NMDS analysis
    Parse labeled and unlabeled data into CSVs that can be used to make a phyloseq object in R.
    Phyloseq object will be used to do NMDS ordinations and PERMANOVAs on CLR-transformed data. 
    
Inputs:
    - PSMs DataFrame for labeled and unlabeled subsets of data
        - Labeled: use Percolator output
        - Unlabeled: use TSV file output by subset_unlabeled.py

Outputs:
    - 2 OTU tables; for both labeled and unlabeld subsets of data
    - 2 Taxonomy tables; for both labeled and unlabeld subsets of data
    - 2 Metadata tables; for both labeled and unlabeld subsets of data

Usage:
    python generate_phyloseq.py \
        -lin [percolator output] \
        -uin [subset_unlabeled.py output] \
        -n [sample names lookup table] \
        -t [MGYG metadata] \
        -ml [Labeled metadata output] \
        -mu [Unlabeled metadata output] \
        -tl [Labeled Taxonomy table output] \
        -tu [Unlbeled Taxonomy table output] \
        -ol [Labeled OTU table output] \
        -ou [Unlabeled OTU table output] \

Notes:
    ASVs: genera names
    Abundance measure: unlabeled or labeled spectral counts
"""

import pandas as pd
import numpy as np
import argparse

class parseMetadata():
    """
    Parse sample naming metadata and MGYG metadata 

    Generate dictionaries to query sample names, label status, genera names, 
    and entire lineage of each genus

    Attributes:
        namesDict (file path) : path to sample naming metadata 
        MGYGMeta (file path) : path to MGYG taxonomy metadata 
    """

    def __init__(self, namesDict, MGYGMeta):
        self.namesDict = namesDict
        self.MGYGMeta = MGYGMeta

    def sampleMetadata(self):
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
        sampleLookup = pd.read_csv(self.namesDict)
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
    
    def parseMGYGData(self):
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
        metadataDf = pd.read_csv(self.MGYGMeta, sep = '\t', header = 0, usecols = ['Genome', 'Lineage'])
        lineageDict = {}
        fullLineageDict = {}
        for isolate, lineage in metadataDf.itertuples(index = False):
            splitLineage = lineage.split(';')
            ### splitLineage is a list that can be indexed to get any taxonomic rank in the lineage
            ### Index 5 corresponds to genera 
            lineageDict[isolate] = splitLineage[5]
            fullLineageDict[isolate] = splitLineage
        return lineageDict, fullLineageDict

class generatePhyloseqObjData():
    """
    Parse data and generate data structures that can be used to make a phyloseq object

    Attributes:
        sampleNamingDict (dict) :  Dictionary to convert MS/MS file name to sample name
        sampleStatusDict (dict) : Dictionary that can be used to get label status of each sample (labeled/unlabeled)
        isolateDict (dict) : Lookup dict to get genus name for each taxon ID
        isolateLineageDict (dict) : Lookup dict to get full taxonomic lineage for each taxon ID
    """

    def __init__(self, sampleNamingDict, sampleStatusDict, isolateDict, isolateLineageDict):
        self.sampleNamingDict = sampleNamingDict
        self.sampleStatusDict = sampleStatusDict
        self.isolateDict = isolateDict
        self.isolateLineageDict = isolateLineageDict

    def parseWithTaxonomy(self, fileName, status):
        """
        Parse unlabeled or labeled PSM data 
        Insert genus names and taxonomic lineages 

        Parameters
        ----------
        fileName : file path
            TSV file with unlabeled or Labeled PSM data 

        status : str
            String to indicate if TSV contains labeled ('L') or unlabeled ('U') data

        Returns 
        -------
        OTUTable : pandas.DataFrame
            OTU table for phyloseq object

        taxonomyTable : pandas.DataFrame
            Taxonomy table for phyloseq object
        """
        df = pd.read_csv(fileName, sep = '\t', usecols = ['PSMId', 'Proteins', 'FileName'])

        ### Need to remove unlabeled samples from data when making labeled OTU and taxonomy tables
        if status == 'L':
            df['Status'] = df['FileName'].map(self.sampleStatusDict)
            df = df[df["Status"] != "Unlabeled"]
            df = df.drop('Status', axis = 1)

        abundData = []
        taxData = []
        for psm, protein, sample in df.itertuples(index = False):
            stripProtein = protein.lstrip('{').rstrip('}')
            sampName = self.sampleNamingDict.get(sample)
            if stripProtein.startswith('MGYG'):
                splitProtein = stripProtein.split(',')[0].split('_')
                taxonName = self.isolateDict.get(splitProtein[0])
                taxonLineage = self.isolateLineageDict.get(splitProtein[0])
                abundData.append([taxonName, sampName])
                taxData.append([taxonName, *taxonLineage])

        otuTablsDf = pd.DataFrame(abundData, columns = ['Taxon', 'Sample'])
        ### Assume each row = 1 spectral count
        gbOTUCountDf = otuTablsDf.groupby(['Sample', 'Taxon'], sort = False).size().reset_index()
        OTUTable = pd.pivot_table(gbOTUCountDf, index = 'Taxon', columns = 'Sample', aggfunc = 'sum', sort = False).fillna(0)     
        OTUTable.columns = OTUTable.columns.droplevel()

        taxonomyTable = pd.DataFrame(taxData, columns = ['otu', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']).sort_values(by = 'otu').drop_duplicates('otu')
        return OTUTable, taxonomyTable

    def genPhyMetadata(self, labelStatus):
        """
        Generate sample metadata for phyloseq object
        Independent variables included: treatment group, sample type (cecum or stool), sampling time point, and replicate 

        Parameters
        ----------
        labelStatus : str
            String to indicate if we're making a metadata table for unlabeled or labeled samples

        Returns
        -------
        metadataDf : pandas.DataFrame
            Metadata table for phyloseq object
        """
        metadata = []
        for fileName, sampleName in self.sampleNamingDict.items():
            splitName = sampleName.split('.')
            group = splitName[0]
            sampleType = splitName[0][0]
            replicate = splitName[1]
            time = splitName[0][1::]
            metadata.append([sampleName, group, sampleType, replicate, time])
        metadataDf = pd.DataFrame(metadata, columns = ['Sample', 'Group', 'SampleType', 'Replicate', 'Time']).drop_duplicates('Sample').set_index('Sample')
        
        ### Can't include unlabeled (0 hour) samples in labeled metadata table
        if labelStatus == 'L':
            metadataDf = metadataDf[metadataDf["Time"] != "0"]
        return metadataDf

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-lin', '--labeledDataFrame')
    parser.add_argument('-uin', '--unlabeledDataFrame')
    parser.add_argument('-n', '--names')
    parser.add_argument('-t', '--MGYGMetadata')

    parser.add_argument('-ml', '--labeledMetadata')
    parser.add_argument('-mu', '--unlabeledMetadata')
    parser.add_argument('-tl', '--labeledTaxonomyTab')
    parser.add_argument('-tu', '--unlabeledTaxonomyTab')
    parser.add_argument('-ol', '--labeledOTUTab')
    parser.add_argument('-ou', '--unlabeledOTUTab')
    args = parser.parse_args()

    metadataParser = parseMetadata(args.names, args.MGYGMetadata)
    sDict, sstatDict = metadataParser.sampleMetadata()
    linDict, flinDict = metadataParser.parseMGYGData()

    genPhyObject = generatePhyloseqObjData(sDict, sstatDict, linDict, flinDict)
    unlab_OTUTAB, unlab_TAXTAB = genPhyObject.parseWithTaxonomy(args.labeledDataFrame, 'L')
    lab_OTUTAB, lab_TAXTAB = genPhyObject.parseWithTaxonomy(args.unlabeledDataFrame, 'U')
    lab_METATAB = genPhyObject.genPhyMetadata('L')
    unlabMETATAB = genPhyObject.genPhyMetadata('U')

    unlab_OTUTAB.to_csv(args.unlabeledOTUTab)
    unlab_TAXTAB.to_csv(args.unlabeledTaxonomyTab, index = False)
    lab_OTUTAB.to_csv(args.labeledOTUTab)
    lab_TAXTAB.to_csv(args.labeledTaxonomyTab, index = False)
    lab_METATAB.to_csv(args.labeledMetadata)
    unlabMETATAB.to_csv(args.unlabeledMetadata)

if __name__ == "__main__":
    main()