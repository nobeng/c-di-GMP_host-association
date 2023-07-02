# Python script to sort and index .bam-files, then deduplicate reads
# !/usr/bin/env python

# Load package
import os
import subprocess

# Define directories and samples
dirAlignments = "/sfs/fs5/home-sh/sunzm524/genomesMorphotypes/genomesMorphotypesAligned"
dirPicardTools = "/sfs/fs5/home-sh/sunzm524/software/picard-2.22.2"

# Collect unique sample names in directory
allSamples = os.listdir('{}'.format(dirAlignments)) # List files in dir.
allSamples = [i.split('_paired.bam', 1)[0] for i in allSamples] # Collect sample names before read identifier

# Define function to sort and index bam files
def sortNindexBAM(sampleName):
    sort_cmd = "samtools sort {}/{}_paired.bam -o {}/{}_paired_sort.bam".format(dirAlignments, sampleName, dirAlignments, sampleName)
    subprocess.call(sort_cmd, shell=True)

    index_cmd = "samtools index {}/{}_paired_sort.bam".format(dirAlignments, sampleName)
    subprocess.call(index_cmd, shell = True)

# Define function to identify and delete duplicate sequences
def markDuplicatesAndRemove(sampleName):
    deduplicate_cmd = "java -jar {}/picard.jar MarkDuplicates INPUT={}/{}_paired_sort.bam OUTPUT={}/{}_deduplicated.bam M={}/{}_marked_dup_metrics.txt REMOVE_DUPLICATES=true".format(dirPicardTools, dirAlignments, sampleName, dirAlignments, sampleName, dirAlignments, sampleName)
    subprocess.call(deduplicate_cmd, shell=True)

# Call functions, i.e. run function fastQC on command line
if __name__=='__main__':
    for s in allSamples:
         sortNindexBAM(s)

    for s in allSamples:
         markDuplicatesAndRemove(s)
    
    	 
