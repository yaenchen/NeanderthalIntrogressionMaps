#!/bin/bash

HOMEDIR='/wynton/home/capra/ychen39/'
BEDTOOLSDIR='/wynton/home/cbi/shared/software/CBI/bedtools2-2.31.1/bin'

## 0. In case TMPDIR is not set, set it to local /scratch, 
scratch_dir=/wynton/scratch/ychen39/CEU_CHBS/

# get variables from system input in bash


module load CBI miniforge3

# load conda environment
conda activate jupyter

${BEDTOOLSDIR}/bedtools shuffle -i ${scratch_dir}/${algorithm_name}_${person} -g ${HOMEDIR}/introgression_methods/data/human_hg19_gaps_removed_common_chr_autosomes_no_prefix.bed > ${scratch_dir}/${algorithm_name}_Shuffled_${person}_${SGE_TASK_ID}.bg
# sort
sort -k1,1 -k2,2n ${scratch_dir}/${algorithm_name}_Shuffled_${person}_${SGE_TASK_ID}.bg -o ${scratch_dir}/${algorithm_name}_Shuffled_${person}_${SGE_TASK_ID}.bg
