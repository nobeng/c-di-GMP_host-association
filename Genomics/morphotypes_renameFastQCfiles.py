##################################################################################
# Python3 code to rename fastqc files based on sample identity and sequencing ID #
#                          Nancy Obeng, Jan. 2020                                #
##################################################################################
# !/usr/bin/env python
# Load packages
import os
import csv
import re

# Note: this script works well with this legend. However, it has magic numbers for the columns to select in the legend.

# Define directories
dirRawSeq="/sfs/fs5/home-sh/sunzm524/genomesMorphotypes"
legendSeqNames = '{}/'.format(dirRawSeq) + "201910_sequencing_morphotypes_pipetOrder.Conc.csv"

# Define functions
def renameUsingLegend():
    new_name = [] # Empty arrays to collect old and new labels
    old_name =[]

    # Collect legend in two arrays with the sequencing IDs and the corresponding and the exp. label, respectively
    with open("{}".format(legendSeqNames)) as csvDataFile:
        csvReader = csv.reader(csvDataFile)  # Open csv-file
        for row in csvReader:
            old_name.append(row[4]) # Collect sample names from the 4th (5th) column w/o the title
            new_name.append(row[5]) # Collect new labels from the 5th (6th) column w/o the title

    # Collect list of files in the directory with the sequencing data
    my_files = list(filter(lambda x:'.fastq.gz' in x, os.listdir("{}".format(dirRawSeq))))

    # Go through files and re-label according to the legend
    for f in my_files:
        # Extract the sequencing ID by splitting the file name
        seqID = re.split('L1_|_', f)[1] # Specifically, split after "L1_" and before "_", then take the second entry, i.e. index 1
        fileEnding = re.split('L001', f)[1] # Collect end of file name including read info

        # Find position of sequencing ID in the legend and find label accordingly
        positionInLegend = old_name.index(seqID)
        newSeqLabel = new_name[positionInLegend]

        dst = newSeqLabel + fileEnding # new name
        src = "{}/".format(dirRawSeq) + f # old name with dir
        dst = "{}/".format(dirRawSeq) + dst # new name with dir

        os.rename(src, dst)

# Driver code
if __name__=='__main__':
    renameUsingLegend();
