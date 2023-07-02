# Python script to conver .bcf to vcf
# !/usr/bin/env python

# Import packages
import subprocess
import os

# Define directories and samples
dirBCFtools = "/sfs/fs5/home-sh/sunzm524/software/bcftools"
dirVarCalls = "/sfs/fs5/home-sh/sunzm524/genomesMorphotypes/genomesMorphotypesVariants"

allSamples = os.listdir('{}'.format(dirVarCalls)) # List files in dir.
allSamples = list(filter(lambda x:'_var_filtered.bcf' in x, allSamples))
allSamples = [i.split('_var', 1)[0] for i in allSamples] # Collect sample names before read identifier
uniqueSamples = list(set(allSamples))

# Explanation for parameters used below
# Define function to call variants
def changeFileFormat(sampleName):
    changeFile_cmd = "{}/bcftools convert {}/{}_var_filtered.bcf -o {}/{}_var_filtered.vcf".format(dirBCFtools, dirVarCalls, sampleName, dirVarCalls, sampleName)
    subprocess.call(changeFile_cmd, shell=True)

# Call functions
if __name__ == '__main__':
    for s in uniqueSamples:
        changeFileFormat(s);