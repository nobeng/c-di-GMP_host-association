#!/bin/bash
#PBS -b 2
#PBS -l cpunum_job=32
#PBS -l elapstim_req=02:00:00
#PBS -l memsz_job=30gb
#PBS -N varCallMorphotypes
#PBS -o varCallMorphotypes.out
#PBS -e varCallMorphotypes.err
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

# Start the computation4
python3 morphotypeCallVar.py

# Output of used resources (computation time, main memory) after the job
/opt/nec/nqsv/bin/qstat -f ${PBS_JOBID/0:/}