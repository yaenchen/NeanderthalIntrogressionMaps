#!/bin/bash
#$ -M yaen.chen@ucsf.edu
#$ -m a
#$ -o /wynton/home/capra/ychen39/introgression_methods/introgression_tools/bin/job_scripts/ceu_individual_jaccard/ceu_individual_jaccard.o
#$ -e /wynton/home/capra/ychen39/introgression_methods/introgression_tools/bin/job_scripts/ceu_individual_jaccard/ceu_individual_jaccard.e
#$ -l h_rt=5:00:00
#$ -l mem_free=50G
#$ -l scratch=100G


#qsub /wynton/home/capra/ychen39/introgression_methods/introgression_tools/bin/job_scripts/ceu_individual_jaccard/ceu_individual_jaccard.sh

HOMEDIR='/wynton/home/capra/ychen39/'
BEDTOOLSDIR='/wynton/home/cbi/shared/software/CBI/bedtools2-2.31.1/bin'

## 0. In case TMPDIR is not set, set it to local /scratch, 
##    if it exists, otherwise to /tmp
if [[ -z "$TMPDIR" ]]; then
  if [[ -d /scratch ]]; then TMPDIR=/scratch/$USER; else TMPDIR=/tmp/$USER; fi
  mkdir -p "$TMPDIR"
  export TMPDIR
fi

echo $TMPDIR

module load CBI miniforge3

# load conda environment
conda activate jupyter

python3 /wynton/home/capra/ychen39/introgression_methods/introgression_tools/bin/job_scripts/ceu_individual_jaccard/ceu_individual_jaccard.py $TMPDIR