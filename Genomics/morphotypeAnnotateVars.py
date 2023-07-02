# Python script to call Variants using samtools and BCFtools
# !/usr/bin/env python

# Import packages
import subprocess
import os

# Define directories
dirVarCalls = "/sfs/fs5/home-sh/sunzm524/genomesMorphotypes/genomesMorphotypesVariants"
refDatabase = "myb11"

# Collect sample names
allSamples = os.listdir('{}'.format(dirVarCalls)) # List files in dir.
allSamples = list(filter(lambda x:'_var_filtered.vcf' in x, allSamples))  # Collect all fastq.gz files
allSamples = [i.split('_var_filtered', 1)[0] for i in allSamples] # Collect sample names before read identifier

# Define function to annotate called variants using SnpEff
def annotateUsingSnpEff(sampleName):
    renameChromosome_cmd = "cat {}/{}_var_filtered.vcf | sed s/^NZ_CP023272.1/NZ_CP023272/ > {}/{}_updated_chr.vcf".format(dirVarCalls, sampleName, dirVarCalls, sampleName)
    subprocess.call(renameChromosome_cmd, shell = True)

    annotateVar_cmd = "java -jar snpEff.jar -geneId {} {}/{}_updated_chr.vcf > {}/{}_var_ann.vcf -stats {}/{}_stats".format(refDatabase, dirVarCalls, sampleName, dirVarCalls, sampleName, dirVarCalls, sampleName)
    subprocess.call(annotateVar_cmd, shell = True)

    filterByQuality_cmd = "java -jar snpSift.jar filter 'QUAL >= 20' {}/{}_var_ann.vcf > {}/{}_var_ann_filtered.vcf".format(dirVarCalls, sampleName, dirVarCalls, sampleName)
    subprocess.call(filterByQuality_cmd, shell = True)

# Call function
if __name__ == '__main__':
    for s in allSamples:
        annotateUsingSnpEff(s)