# Python script to call Variants using VarScan2
# !/usr/bin/env python

# Import packages
import subprocess
import os

# Define directories and samples
dirAlignments = "/sfs/fs5/home-sh/sunzm524/genomesMorphotypes/genomesMorphotypesAligned"
dirVarCalls = "/sfs/fs5/home-sh/sunzm524/genomesMorphotypes/genomesMorphotypesVariants"
dirVarScan = "/sfs/fs5/home-sh/sunzm524/software"

dirReference = "/sfs/fs5/home-sh/sunzm524/ancestralGenomes"
dirPloidyFile = dirReference
nameReference = "myb11_ancestralREFSEQ.fasta"

allSamples = os.listdir('{}'.format(dirAlignments)) # List files in dir.
allSamples = list(filter(lambda x:'_deduplicated.bam' in x, allSamples))  # Collect all fastq.gz files
allSamples = [i.split('_deduplicated', 1)[0] for i in allSamples] # Collect sample names before read identifier
uniqueSamples = list(set(allSamples))
ancestralSample = "MT48_S16"
uniqueSamples.remove(ancestralSample) # given comparison with ancestor, remove it from list of evolved samples

# Explanation for parameters used below
# Samtools
# -u --> continue with uncompressed file
# -f --> input fasta file

# Function to collect copy number
def callCopyNumber(sampleName):
    # Collect copy numbers
    callCopyNumberStep1_cmd = "samtools mpileup -B -q 1 -f   {}/{} {}/{}_deduplicated.bam {}/{}_paired_sort.bam | java -jar {}/VarScan.v2.3.9.jar copynumber varScan {}/{} --mpileup 1 ".format(dirReference, nameReference, dirAlignments, ancestralSample, dirAlignments, sampleName, dirVarScan, dirVarCalls, sampleName)
    subprocess.call(callCopyNumberStep1_cmd, shell=True)

    # Correct for GC content
    callCopyNumberStep2_cmd = "java -jar {}/VarScan.v2.3.9.jar copyCaller {}/{}.copynumber --output-file {}/{}_copyNumsCalled".format(dirVarScan, dirVarCalls, sampleName, dirVarCalls, sampleName)
    subprocess.call(callCopyNumberStep2_cmd, shell=True)

# Call functions
if __name__ == '__main__':
    for s in uniqueSamples:
        callCopyNumber(s);