#!/bin/bash
#$ -M yaen.chen@ucsf.edu
#$ -m a
#$ -o /wynton/home/capra/ychen39/introgression_methods/bin/output/bootstrap_bstat_dist.o
#$ -e /wynton/home/capra/ychen39/introgression_methods/bin/output/bootstrap_bstat_dist.e
#$ -l h_rt=48:00:00
#$ -l mem_free=50G


#qsub /wynton/home/capra/ychen39/introgression_methods/bin/bootstrap_bstat_dist.sh

HOMEDIR='/wynton/home/capra/ychen39/'
BEDTOOLSDIR='/wynton/home/cbi/shared/software/CBI/bedtools2-2.31.1/bin'

module load CBI miniforge3

# load conda environment
conda activate jupyter

python3 /wynton/home/capra/ychen39/introgression_methods/bin/bootstrap_bstat_dist.py