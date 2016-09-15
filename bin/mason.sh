#!/bin/bash
#PBS -k o
#PBS -l nodes=1:ppn=8,vmem=100gb,walltime=72:00:00
#PBS -M wrshoema@indiana.edu
#PBS -m abe
#PBS -j oe
module load python

python /N/dc2/projects/Lennon_Sequences/2016_MutationDorm/Python/test_mason.py
