# Python3 script to map paired, trimmed morphotype sequences to the ancestral MYb11 genome
# !/usr/bin/env python

# General notes on settings:
# - end-to-end mode is used (no clipping of sequences)
#   - using very-sensitive settings (using the pre-set alignment options)
# - default: upstream/downstream mate orientations for a valid paired-end alignment against the forward reference strand
# - input reads are fastq files (--> -q)
# - output files are .sam-files (--> -S)
# - input qualities are ASCII characters following the current Illumina framework (--> --phred33)
# - CPUs to use: 8 (--> -p)

# Load packages
import os
import subprocess

# Define directories and variables
dirTrimmedSeq="/sfs/fs5/home-sh/sunzm524/genomesMorphotypes/genomesMorphotypesTrimmed"
dirIndexedRef = "/sfs/fs5/home-sh/sunzm524/ancestralGenomes"
nameIndexedRef = "index_myb11"
dirAlignments = "/sfs/fs5/home-sh/sunzm524/genomesMorphotypes/genomesMorphotypesAligned"

# Collect unique sample names in directory
allSamples = os.listdir('{}'.format(dirTrimmedSeq)) # List files in dir.
allSamples = list(filter(lambda x:'.fq.gz' in x, allSamples))  # Collect all fastq.gz files
allSamples = [i.split('_R', 1)[0] for i in allSamples] # Collect sample names before read identifier
uniqueSamples = list(set(allSamples))

sampleName = "MT48_S16"

# Define function to map reads using bowtie2
def mapUsingbowtie2(sampleName):
    fileName = sampleName
    forwardPair = "{}/{}_R1_paired.fq.gz".format(dirTrimmedSeq, sampleName)
    reversePair = "{}/{}_R2_paired.fq.gz".format(dirTrimmedSeq, sampleName)

    run_bowtie2_cmd = "bowtie2 -p 8 --very-sensitive -q --phred33 -x {}/{} -1 {} -2 {} |samtools view -bS -> {}/{}_paired.bam ".format(dirIndexedRef, nameIndexedRef, forwardPair, reversePair, dirAlignments, fileName)
    subprocess.call(run_bowtie2_cmd, shell=True)

# Driver code
if __name__=='__main__':
    for s in uniqueSamples:
        mapUsingbowtie2(s)
