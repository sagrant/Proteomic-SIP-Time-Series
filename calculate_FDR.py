#!/usr/bin/env python3
import pandas as pd
import numpy as np
import argparse

'''
Calculate FDR for entire dataset
FDR = ((total labeled PSMs in control samples / total spectra in control samples) * total spectra in labeled samples / labeled PSMs in labeled samples)

call with: 
python calculate_FDR.py -n [sample name dict] -p [Percolator output] -s [CSV with counts of spectra aquired in each sample]

'''
class calculateDataSetFDR():

    def __init__(self, SIPPsms, sampleSpectra):
        self.SIPPsms = SIPPsms
        self.sampleSpectra = sampleSpectra

    def sampleMetadata(self, namesDictIn):
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
        gbSIPdf = self.SIPPsms.groupby(self.SIPPsms.iloc[:,2])
        totLabeled_unlabeledSamples = []
        totLabeled_labeledSamples = []
        for x, sampleGroup in gbSIPdf:
            sampleName = sDict.get(x)
            status = ssDict.get(sampleName)
            if status == 'Unlabeled':
                totLabeled_unlabeledSamples.append(len(sampleGroup))
            else:
                totLabeled_labeledSamples.append(len(sampleGroup))
        return totLabeled_unlabeledSamples, totLabeled_labeledSamples

    def SpectraCounts(self, ssDict):
        self.sampleSpectra['Status'] = self.sampleSpectra.iloc[:,0].map(ssDict)
        unlabeledSampleSpectra = self.sampleSpectra.set_index('Status').loc['Unlabeled']
        totSpectra_unlabeledSamples = unlabeledSampleSpectra.iloc[:,1].sum()
        labeledSampleSpectra = self.sampleSpectra.set_index('Status').loc['Labeled']
        totSpectra_labeledSamples = labeledSampleSpectra.iloc[:,1].sum()
        return totSpectra_unlabeledSamples, totSpectra_labeledSamples

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--namesDict')
    parser.add_argument('-p', '--percolatorOutput')
    parser.add_argument('-s', '--spectra')
    args = parser.parse_args()

    SIPdf = pd.read_csv(args.percolatorOutput, sep = '\t', header = 0, usecols = [0, 25, 26])
    spectraDf = pd.read_csv(args.spectra).dropna()

    calcFDR = calculateDataSetFDR(SIPdf, spectraDf)
    sampleDict, statusDict = calcFDR.sampleMetadata(args.namesDict)
    totalLabeled_unlabeledSamples, totalLabeled_labeledSamples = calcFDR.PSMCounts(sampleDict, statusDict)
    totalSpectra_unlabeledSamples, totalSpectra_labeledSamples = calcFDR.SpectraCounts(statusDict)

    FDR = ((sum(totalLabeled_unlabeledSamples) / totalSpectra_unlabeledSamples) * totalSpectra_labeledSamples) / sum(totalLabeled_labeledSamples)
    print('------------------------------------------')
    print(f'FDR for this dataset = {FDR}')
    print('------------------------------------------')
    print('------------------------------------------')

if __name__ == "__main__":
    main()