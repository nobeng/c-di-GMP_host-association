#!/bin/bash
#PBS -b 2
#PBS -l cpunum_job=32
#PBS -l elapstim_req=02:00:00
#PBS -l memsz_job=15gb
#PBS -N myb11_copyNum
#PBS -o myb11_copyNum.out
#PBS -j o
#PBS -q clexpress
#PBS -l cputim_job=32:00:00
#PBS -m abe
#PBS -M nobeng@zoologie.uni-kiel.de

# Initialize Python3
module load python3.7.4
module load samtools1.5
module load java1.8.0

# Change into subdirectory containing the Python script
cd $HOME/scripts/pythonScripts

# Start the computation
python3 morphotypes_copyNumber.py && python3 popsMYb11_copyNumber.py

# Output of used resources (computation time, main memory) after the job
/opt/nec/nqsv/bin/qstat -f ${PBS_JOBID/0:/}