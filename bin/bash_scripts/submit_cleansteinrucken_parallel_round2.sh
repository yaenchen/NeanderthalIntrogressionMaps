#!/bin/bash
#$ -M yaen.chen@ucsf.edu
#$ -m a
#$ -o /wynton/home/capra/ychen39/introgression_methods/bin/output/submit_cleansteinrucken_parallel_round2.o
#$ -e /wynton/home/capra/ychen39/introgression_methods/bin/output/submit_cleansteinrucken_parallel_round2.e

# go to directory containing files
cd /wynton/group/capra/data/introgression_maps/steinruecken_2018/CalledTracts_March2018/
# get list of tar.gz files into targz_files variable
targz_files=`ls *.tar.gz`

# make round 2 output directory
mkdir /wynton/home/capra/ychen39/introgression_methods/cleaned/steinruecken_2018/round2/

for targz_file in $targz_files
do
    qsub -v targz=$targz_file /wynton/home/capra/ychen39/introgression_methods/bin/run_cleansteinrucken_parallel_round2.sh
done