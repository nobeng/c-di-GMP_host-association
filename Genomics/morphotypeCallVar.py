# Python script to call Variants using samtools and BCFtools
# !/usr/bin/env python

# Import packages
import subprocess
import os

# Define directories and samples
dirAlignments = "/sfs/fs5/home-sh/sunzm524/genomesMorphotypes/genomesMorphotypesAligned"
dirBCFtools = "/sfs/fs5/home-sh/sunzm524/software/bcftools"
dirVarCalls = "/sfs/fs5/home-sh/sunzm524/genomesMorphotypes/genomesMorphotypesVariants"

dirReference = "/sfs/fs5/home-sh/sunzm524/ancestralGenomes"
dirPloidyFile = dirReference
nameReference = "myb11_ancestralREFSEQ.fasta"

allSamples = os.listdir('{}'.format(dirAlignments)) # List files in dir.
allSamples = list(filter(lambda x:'_deduplicated.bam' in x, allSamples))  # Collect all fastq.gz files
allSamples = [i.split('_deduplicated', 1)[0] for i in allSamples] # Collect sample names before read identifier
uniqueSamples = list(set(allSamples))

# Explanation for parameters used below
# -u --> continue with uncompressed file
# -f --> input fasta file

# -c --> consensus caller
# -v --> variants only (don't print matches between ref. and seq.)

# -s --> soft filter with the name LowQual
# -e --> exclude sites for which the Base Alignment Quality is lower than in this case 20

# Define function to call variants
def callVars(sampleName):
    callVars_cmd = "samtools mpileup -u -f {}/{} {}/{}_deduplicated.bam | " \
                  "{}/bcftools call -cv --ploidy-file {}/haploid.ploidy.txt > {}/{}_var_raw.bcf".format(dirReference, nameReference, dirAlignments,sampleName, dirBCFtools,dirPloidyFile, dirVarCalls, sampleName, dirBCFtools,dirVarCalls, sampleName)
    subprocess.call(callVars_cmd, shell=True)

    filterVars_cmd = "{}/bcftools filter {}/{}_var_raw.bcf -s LowQual -e'%QUAL<20 || DP>100' > {}/{}_var_filtered.bcf ".format(dirBCFtools, dirVarCalls, sampleName, dirVarCalls, sampleName)
    subprocess.call(filterVars_cmd, shell=True)

# Call functions
if __name__ == '__main__':
    for s in uniqueSamples:
        callVars(s);