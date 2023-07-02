#!/bin/bash
#PBS -b 2
#PBS -l cpunum_job=32
#PBS -l elapstim_req=02:00:00
#PBS -l memsz_job=15gb
#PBS -N morphotypeQC_trimmed
#PBS -o morphotypeQC_trimmed.out
#PBS -j o
#PBS -q clexpress
#PBS -l cputim_job=32:00:00
#PBS -m abe
#PBS -M nobeng@zoologie.uni-kiel.de

# Initialize Python3
module load python3.7.4

# Change into subdirectory containing the Python script
cd /sfs/fs5/home-sh/sunzm524/scripts/pythonScripts

# Start the computation
python3 qualitycontrol_FastQC_trimmed.py

# Output of used resources (computation time, main memory) after the job
/opt/nec/nqsv/bin/qstat -f ${PBS_JOBID/0:/}