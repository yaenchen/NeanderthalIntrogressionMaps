#!/bin/bash
#$ -M yaen.chen@ucsf.edu
#$ -m a
#$ -o /wynton/home/capra/ychen39/introgression_methods/bin/output/run_cleansteinrucken_parallel_output.o
#$ -e /wynton/home/capra/ychen39/introgression_methods/bin/output/run_cleansteinrucken_parallel_error.e
#$ -l h_rt=48:00:00
#$ -l mem_free=50G

module load CBI miniconda3

# load conda environment
conda activate jupyter

# Run script to generate a membership matrix of introgressed variants detected by various methods
python3 /wynton/home/capra/ychen39/introgression_methods/bin/clean_steinrucken_parallel_round2.py $targz