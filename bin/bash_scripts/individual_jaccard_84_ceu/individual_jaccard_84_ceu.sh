#!/bin/bash
#$ -o /wynton/home/capra/ychen39/introgression_methods/introgression_tools/bin/job_scripts/individual_jaccard_84_ceu/individual_jaccard_84_ceu.o
#$ -e /wynton/home/capra/ychen39/introgression_methods/introgression_tools/bin/job_scripts/individual_jaccard_84_ceu/individual_jaccard_84_ceu.e
#$ -l h_rt=5:00:00
#$ -l mem_free=50G
#$ -l scratch=100G

# qsub /wynton/home/capra/ychen39/introgression_methods/introgression_tools/bin/job_scripts/individual_jaccard_84_ceu/individual_jaccard_84_ceu.sh
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

python3 /wynton/home/capra/ychen39/introgression_methods/introgression_tools/bin/job_scripts/individual_jaccard_84_ceu/individual_jaccard_84_ceu.py $TMPDIR 100 /wynton/home/capra/ychen39/introgression_methods/cleaned/introgression_tools/individual_comparison/CEU84/
# tmp dir
# number of shuffles
# seed for shuffle
# outputdir