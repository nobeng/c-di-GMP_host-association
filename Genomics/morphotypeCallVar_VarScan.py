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

# Explanation for parameters used below
# Samtools
# -u --> continue with uncompressed file
# -f --> input fasta file

# VarScan
# --min-coverage	Minimum read depth at a position to make a call
minCov=10 # based on comminucation with Leifs
# 	--min-reads2	Minimum supporting reads at a position to call variants [2]
minReads=2 # default setting
# 	--min-avg-qual	Minimum base quality at a position to count a read [15]
minAvgQual=20
# 	--min-var-freq	Minimum variant allele frequency threshold [0.01]
minVarFreq=0.1 # communication with Hinrich
# 	--p-value	Default p-value threshold for calling variants [99e-02]
# 	--strand-filter	Ignore variants with >90% support on one strand [1]
# 	--output-vcf	If set to 1, outputs in VCF format
# 	--variants	Report only variant (SNP/indel) positions (mpileup2cns only) [0]

# Define function to call variants
def callVars(sampleName):
    # Call SNPs (VarScan)
    callSNPs_cmd = "samtools mpileup -f {}/{} {}/{}_deduplicated.bam | java -jar {}/VarScan.v2.3.9.jar mpileup2snp --min-coverage {} --min-reads2 {} --min-avg-qual {} --min-var-freq {} --output-vcf 1 > {}/{}_VarScan_SNPs.vcf".format(dirReference, nameReference, dirAlignments, sampleName, dirVarScan, minCov, minReads, minAvgQual, minVarFreq,dirVarCalls, sampleName)
    subprocess.call(callSNPs_cmd, shell=True)

    # Call Indels
    callIndels_cmd = "samtools mpileup -f {}/{} {}/{}_deduplicated.bam | java -jar {}/VarScan.v2.3.9.jar mpileup2indel --mpileup 1 --min-coverage {} --min-reads2 {} --min-avg-qual {} --min-var-freq {} --output-vcf 1 > {}/{}_VarScan_Indels.vcf".format(dirReference, nameReference, dirAlignments, sampleName, dirVarScan, minCov, minReads, minAvgQual, minVarFreq, dirVarCalls, sampleName)
    subprocess.call(callIndels_cmd, shell=True)

# Call functions
if __name__ == '__main__':
    for s in uniqueSamples:    
        callVars(s);