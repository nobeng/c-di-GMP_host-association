#!/bin/bash
#PBS -b 2
#PBS -l cpunum_job=32
#PBS -l elapstim_req=02:00:00
#PBS -l memsz_job=30gb
#PBS -N snpEffMorphotypes
#PBS -o snpEffMorphotypes.out
#PBS -e snpEffMorphotypes.err
#PBS -j o
#PBS -q clexpress
#PBS -l cputim_job=32:00:00
#PBS -m b
#PBS -M nobeng@zoologie.uni-kiel.de

# Initialize Java
module load java1.8.0
module load python3.7.4

# Change into subdirectory containing the called variants
dirSnpEff=/sfs/fs5/home-sh/sunzm524/software/snpEff/
dirPythonScripts=/sfs/fs5/home-sh/sunzm524/scripts/pythonScripts

# Set directory
cd ${dirSnpEff}

# Start the computation
#python3 ${dirPythonScripts}/morphotypeAnnotateVars.py
python3 ${dirPythonScripts}/morphotypeAnnotateVars_VarScan.py

# Output of used resources (computation time, main memory) after the job
/opt/nec/nqsv/bin/qstat -f ${PBS_JOBID/0:/}


