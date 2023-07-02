# Python script to trim paired-end Illumina reads of morphotype samples (in Python3)
# Important: this assumes TruSeq-sequence adapters
# !/usr/bin/env python

# Import packages
import subprocess
import os

# Define directories
dirTrim = "/sfs/fs5/home-sh/sunzm524/software/Trimmomatic-0.39/"
dirRawSeq="/sfs/fs5/home-sh/sunzm524/genomesMorphotypes"
dirTrimmedSeq="/sfs/fs5/home-sh/sunzm524/genomesMorphotypes/genomesMorphotypesTrimmed"
dirFastQC="/sfs/fs5/home-sh/sunzm524/software/FastQC/"
dirOutputFastQC_trim="/sfs/fs5/home-sh/sunzm524/genomesMorphotypes/morphotypeFastQC_trimmed"

# Define input data
allSamples = os.listdir('{}'.format(dirRawSeq)) # List files in dir.
allSamples = list(filter(lambda x:'.fastq.gz' in x, allSamples))  # Collect all fastq.gz files
allSamples = [i.split('_R', 1)[0] for i in allSamples] # Collect sample names before read identifier
uniqueSamples = list(set(allSamples))

# Define global trimmomatic parameters
trimVersion = "trimmomatic-0.39"
phredThreshold = "phred33"
adaptersUsed = "/sfs/fs5/home-sh/sunzm524/software/Trimmomatic-0.39/adapters/NexteraPE-PE.fa" # IlluminaClip
seedMismatches = "2" # commonly used
palindromeClipThreshold = "30" # Usadella lab info
simpleClipThreshold = "10" # interm. val. Usadella lab; using sees mismatch, so no palindrom parameters needed
windowSize = "4" # Usadella lab example
requiredQuality = "20" # Based on FastQC ok quality threshold, more conservative than suggested by Usadella lab
qualityLead = "5" # Leif
qualityTrail = "5" # Leif
lengthHeadcrop = "5" # Leif
lengthMin = "36" # Usadella lab

# Define functions
def runTrimmomatic(sampleName):
    # Define variables: Sequences in and out
    inForward = "{}/{}_R1_001.fastq.gz".format(dirRawSeq, sampleName)  # Paired-end fastq input and reference genome
    inReverse = "{}/{}_R2_001.fastq.gz".format(dirRawSeq, sampleName)  # Paired-end fastq input and reference genome
    outForwardPaired = "{}/{}_R1_paired.fq.gz".format(dirTrimmedSeq, sampleName)
    outForwardUnpaired = "{}/{}_R1_unpaired.fq.gz".format(dirTrimmedSeq, sampleName)
    outReversePaired = "{}/{}_R2_paired.fq.gz".format(dirTrimmedSeq, sampleName)
    outReverseUnpaired = "{}/{}_R2_unpaired.fq.gz".format(dirTrimmedSeq, sampleName)

    # Run function with given parameters
    trim_cmd = "java -jar {}/{}.jar PE -{} {} {} {} {} {} {} ILLUMINACLIP:{}:{}:{}:{} SLIDINGWINDOW:{}:{} LEADING:{} TRAILING:{} HEADCROP:{} MINLEN:{}".format(dirTrim, trimVersion, phredThreshold, inForward, inReverse, outForwardPaired, outForwardUnpaired, outReversePaired, outReverseUnpaired, adaptersUsed, seedMismatches, palindromeClipThreshold, simpleClipThreshold, windowSize, requiredQuality, qualityLead, qualityTrail, lengthHeadcrop, lengthMin)
    subprocess.call(trim_cmd, shell=True)


# Call functions, i.e. run function fastQC on command line
if __name__=='__main__':
    # Trim all available sequences
    for s in uniqueSamples:
        runTrimmomatic(s)

    # Find all trimmed sequences
    allSamples = os.listdir('{}'.format(dirTrimmedSeq))  # List files in dir.
    unpairedSeqs = list(filter(lambda x: 'unpaired.fastq.gz' in x, allSamples))  # Collect all unpaired sequences

    # Remove all files that contain unpaired sequences
    for s in unpairedSeqs:
        os.remove("{}/{}".format(dirTrimmedSeq, s))
