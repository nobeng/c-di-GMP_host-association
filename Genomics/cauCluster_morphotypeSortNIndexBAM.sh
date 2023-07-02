#!/bin/bash
#PBS -b 2
#PBS -l cpunum_job=32
#PBS -l elapstim_req=02:00:00
#PBS -l memsz_job=150gb
#PBS -N morphotypes_sortNIndex
#PBS -o morphotypes_sortNIndex.out
#PBS -e morphotypes_sortNIndex.err
#PBS -j o
#PBS -q clexpress
#PBS -l cputim_job=32:00:00
#PBS -m b
#PBS -M nobeng@zoologie.uni-kiel.de

# Initialize Python3
module load python3.7.4
module load samtools1.5

# Change into subdirectory containing the Python script
cd /sfs/fs5/home-sh/sunzm524/scripts/pythonScripts

# Start the computation
python3 morphotype_deduplicate_sortNindexBAM.py

# Output of used resources (computation time, main memory) after the job
/opt/nec/nqsv/bin/qstat -f ${PBS_JOBID/0:/}