#!/bin/bash

# set up directories
HOMEDIR=/wynton/home/capra/ychen39/
INPUTDIR=${HOMEDIR}/introgression_methods/data/vernot_2016/introgressed_tag_snp_frequencies
OUTPUTDIR=${HOMEDIR}/introgression_methods/cleaned/introgression_tools/vernot_2016/

# for SNPs files (not extended LD)
# add data from each population file
for i in ASN EUR PNG SAS; do
    awk -v pop="$i" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $14 "\t" $15 "\t" $16 "\t" pop}' \
    $INPUTDIR/all_tag_snps."$i".merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.bed >> $OUTPUTDIR/neanderthal_introgressed_variants_pops.bed # correct number of snps, should be 526,225
done
# sort by chromosome
sort -u $OUTPUTDIR/neanderthal_introgressed_variants_pops.bed | sort -k 1,1 -k2,2n > $OUTPUTDIR/tmp; mv $OUTPUTDIR/tmp $OUTPUTDIR/neanderthal_introgressed_variants_pops.bed
# add column names
sed -i '1s/^/Chromosome\tStart\tStop\tAncestralAllele\tDerivedAllele\tNeanderthalBase\tDenisovaBase\tHaplotypeTag\tAncestry\n/' $OUTPUTDIR/neanderthal_introgressed_variants_pops.bed

# for extended LD files
# add data from each population file
for i in ASN EUR PNG SAS; do
    awk -v pop="$i" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" pop}' \
    $INPUTDIR/all_tag_snps."$i".merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.bed.extended_LD.sorted >> $OUTPUTDIR/neanderthal_introgressed_variants_pops_extendedLD.bed
done
# sort by chromosome
sort -u $OUTPUTDIR/neanderthal_introgressed_variants_pops_extendedLD.bed | sort -k 1,1 -k2,2n > $OUTPUTDIR/tmp; mv $OUTPUTDIR/tmp $OUTPUTDIR/neanderthal_introgressed_variants_pops_extendedLD.bed
# add column names
sed -i '1s/^/Chromosome\tStart\tStop\tHaplotypeTag\tAncestry\n/' $OUTPUTDIR/neanderthal_introgressed_variants_pops_extendedLD.bed


# for median SNPs + extended LD
# add data from each population file
for i in ASN EUR PNG SAS; do
    awk -v pop="$i" '{print $1 "\t" $2 "\t" $3 "\t" $11 "\t" $12 "\t" pop}' \
    $INPUTDIR/all_tag_snps."$i".merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.median_af.bed.extended_LD >> $OUTPUTDIR/neanderthal_introgressed_variants_pops_medianextendedLD.bed
done
# sort by chromosome
sort -u $OUTPUTDIR/neanderthal_introgressed_variants_pops_medianextendedLD.bed | sort -k 1,1 -k2,2n > $OUTPUTDIR/tmp; mv $OUTPUTDIR/tmp $OUTPUTDIR/neanderthal_introgressed_variants_pops_medianextendedLD.bed
# add column names
sed -i 'Chromosome\tStart\tStop\tSNPsOnHaplotype\tLength\tAncestry\n/' $OUTPUTDIR/neanderthal_introgressed_variants_pops_medianextendedLD.bed
