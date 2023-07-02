#!/bin/bash
#PBS -b 2
#PBS -l cpunum_job=32
#PBS -l elapstim_req=02:00:00
#PBS -l memsz_job=15gb
#PBS -N trimMorphotypes
#PBS -o trimMorphotypes.out
#PBS -j o
#PBS -q clexpress
#PBS -l cputim_job=32:00:00
#PBS -m abe
#PBS -M nobeng@zoologie.uni-kiel.de

# Initialize Python3
module load java1.8.0
module load python3.7.4

# Change into subdirectory containing the Python script
cd /sfs/fs5/home-sh/sunzm524/scripts/pythonScripts

# Start the computation
python3 morphotypeTrimming.py

# Output of used resources (computation time, main memory) after the job
/opt/nec/nqsv/bin/qstat -f ${PBS_JOBID/0:/}