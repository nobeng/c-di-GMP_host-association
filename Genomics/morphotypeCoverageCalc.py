# Python script to calculate coverage using Bedtools
# !/usr/bin/env python

# Import packages
import subprocess
import os

# Define directories and samples
dirAlignments = "/sfs/fs5/home-sh/sunzm524/genomesMorphotypes/genomesMorphotypesAligned"
dirBEDtools = "/sfs/fs5/home-sh/sunzm524/software/bedtools2/bin"

allSamples = os.listdir('{}'.format(dirAlignments)) # List files in dir.
allSamples = list(filter(lambda x:'_deduplicated.bam' in x, allSamples))
allSamples = [i.split('_deduplicated.bam', 1)[0] for i in allSamples] # Collect sample names before read identifier
uniqueSamples = list(set(allSamples))

# Define function to calculate coverage
def calcCov(sampleName):
    # Sort .bam files and use as input to bedtools to get genome coverage
    calcCov_cmd = "{}/bedtools genomecov -ibam {}/{}_deduplicated.bam > {}/{}_coverage.txt". format(dirBEDtools, dirAlignments, sampleName, dirAlignments, sampleName)
    subprocess.call(calcCov_cmd, shell = True)

# Call functions
if __name__=='__main__':
    for s in uniqueSamples:
        calcCov(s);