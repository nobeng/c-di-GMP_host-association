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
allSamples = list(filter(lambda x:'_VarScan_SNPs.vcf' in x, allSamples))  # Collect all fastq.gz files
allSamples = [i.split('_VarScan_SNPs.vcf', 1)[0] for i in allSamples] # Collect sample names before read identifier
variantTypes = ["SNPs", "Indels"]

# Define function to annotate called variants using SnpEff
def annotateUsingSnpEff(sampleName, variantType):
    renameChromosome_cmd = "cat {}/{}_VarScan_{}.vcf | sed s/^NZ_CP023272.1/NZ_CP023272/ > {}/{}_VarScan_{}_upd.vcf".format(dirVarCalls, sampleName, variantType, dirVarCalls, sampleName, variantType)
    subprocess.call(renameChromosome_cmd, shell = True)

    annotateVar_cmd = "java -jar snpEff.jar {} {}/{}_VarScan_{}_upd.vcf > {}/{}_var_ann_{}.vcf -stats {}/{}_stats".format(refDatabase, dirVarCalls, sampleName, variantType, dirVarCalls, sampleName, variantType, dirVarCalls, sampleName)
    subprocess.call(annotateVar_cmd, shell = True)

    #filterByQuality_cmd = "java -jar snpSift.jar filter 'QUAL >= 20' {}/{}_var_ann.vcf > {}/{}_var_ann_filtered.vcf".format(dirVarCalls, sampleName, dirVarCalls, sampleName)
    #subprocess.call(filterByQuality_cmd, shell = True)

# Call function
if __name__ == '__main__':
    for s in allSamples:
        for t in variantTypes:
            annotateUsingSnpEff(s, t);