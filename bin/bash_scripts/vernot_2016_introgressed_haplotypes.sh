#!/bin/bash

# set up directories
HOMEDIR=/wynton/home/capra/ychen39/
INPUTDIR=${HOMEDIR}/introgression_methods/data/vernot_2016/introgressed_haplotypes
OUTPUTDIR=${HOMEDIR}/introgression_methods/cleaned/introgression_tools/vernot_2016
BEDTOOLSDIR=/wynton/home/cbi/shared/software/CBI/bedtools2-2.31.1/bin

# for neanderthal
for i in EAS EUR PNG SAS; do
    # add data from each population file
    awk -v pop="$i" '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" pop}' \
    $INPUTDIR/LL.callset"$i".mr_0.99.neand_calls_by_hap.bed.merged.by_chr.bed >> $OUTPUTDIR/neanderthal_introgressed_haplotypes_individual_no_scores.bed
    # also take the raw data (with scores) and convert to bedgraph so we can use unionbedg
    awk -v pop="$i" '{print $4 "\t" $2 "\t" $3 "\t" $10}' \
    $INPUTDIR/LL.callset"$i".mr_0.99.neand_calls_by_hap.bed >> $OUTPUTDIR/neanderthal_introgressed_haplotypes_individual_scores.bed
    # remove header
    sed -i '1d' $OUTPUTDIR/neanderthal_introgressed_haplotypes_individual_scores.bed
    # run unionbedg
    $BEDTOOLSDIR/bedtools unionbedg -i $OUTPUTDIR/neanderthal_introgressed_haplotypes_individual_no_scores.bed $OUTPUTDIR/neanderthal_introgressed_haplotypes_individual_scores.bed  > $OUTPUTDIR/neanderthal_introgressed_haplotypes_individual.bed

done

# for denisovan
# add data from each population file
for i in PNG; do
    awk -v pop="$i" '{print $1 "\t" $2 "\t" $3 "\t" $5 "\t" pop}' \
    $INPUTDIR/LL.callset"$i".mr_0.99.den_calls_by_hap.bed.merged.by_chr.bed >> $OUTPUTDIR/denisovan_introgressed_haplotypes_individual.bed # correct number of snps, should be 526,225
done

# sort by chromosome
sort -u $OUTPUTDIR/neanderthal_introgressed_haplotypes_individual.bed | sort -k 1,1 -k2,2n > $OUTPUTDIR/tmp; mv $OUTPUTDIR/tmp $OUTPUTDIR/neanderthal_introgressed_haplotypes_individual.bed
sort -u $OUTPUTDIR/denisovan_introgressed_haplotypes_individual.bed | sort -k 1,1 -k2,2n > $OUTPUTDIR/tmp; mv $OUTPUTDIR/tmp $OUTPUTDIR/denisovan_introgressed_haplotypes_individual.bed
