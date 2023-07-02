#!/bin/bash
#PBS -b 2
#PBS -l cpunum_job=32
#PBS -l elapstim_req=02:00:00
#PBS -l memsz_job=15gb
#PBS -N mapMorphotypes_repMT36
#PBS -o mapMorphotypes_repMT36.out
#PBS -e mapMorphotypes_repMT36.err
#PBS -j o
#PBS -q clexpress
#PBS -l cputim_job=32:00:00
#PBS -m b
#PBS -M nobeng@zoologie.uni-kiel.de

# Initialize Python3
module load python3.7.4
module load bowtie2-2.3.3
module load samtools1.5

# Change into subdirectory containing the Python script
cd /sfs/fs5/home-sh/sunzm524/scripts/pythonScripts

# Start the computation
python3 morphotypeMapping_repeat.py

# Output of used resources (computation time, main memory) after the job
/opt/nec/nqsv/bin/qstat -f ${PBS_JOBID/0:/}