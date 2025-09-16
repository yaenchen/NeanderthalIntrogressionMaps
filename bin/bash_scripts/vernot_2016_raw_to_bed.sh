%%bash

# for SNPs files (not extended LD)
# add column names
echo -e "Chromosome\tStart\tStop\tAncestralAllele\tDerivedAllele\tNeanderthalBase\tDenisovaBase\tHaplotypeTag\tAncestry" > $OUTPUTDIR/ancestral_introgressed_variants.bed
# add data from each population file
for i in ASN EUR PNG SAS; do
    awk -v pop="$i" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" $5 "\t" $14 "\t" $15 "\t" $16 "\t" pop}' \
    $INPUTDIR/all_tag_snps."$i".merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.bed >> $OUTPUTDIR/ancestral_introgressed_variants.bed
done
# sort by chromosome
sort -u $OUTPUTDIR/ancestral_introgressed_variants.bed | sort -k 1,1 -k2,2n > $OUTPUTDIR/tmp; mv $OUTPUTDIR/tmp $OUTPUTDIR/ancestral_introgressed_variants.bed

# for extended LD files
# add column names
echo -e "Chromosome\tStart\tStop\tHaplotypeTag\tAncestry" > $OUTPUTDIR/ancestral_introgressed_variants_extendedLD.bed
# add data from each population file
for i in ASN EUR PNG SAS; do
    awk -v pop="$i" '{print $1 "\t" $2 "\t" $3 "\t" $4 "\t" pop}' \
    $INPUTDIR/all_tag_snps."$i".merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.bed.extended_LD.sorted >> $OUTPUTDIR/ancestral_introgressed_variants_extendedLD.bed
done
# sort by chromosome
sort -u $OUTPUTDIR/ancestral_introgressed_variants_extendedLD.bed | sort -k 1,1 -k2,2n > $OUTPUTDIR/tmp; mv $OUTPUTDIR/tmp $OUTPUTDIR/ancestral_introgressed_variants_extendedLD.bed

# for median SNPs + extended LD
# add column names
echo -e "Chromosome\tStart\tStop\tSNPsOnHaplotype\tLength\tAncestry" > $OUTPUTDIR/ancestral_introgressed_variants_medianextendedLD.bed
# add data from each population file
for i in ASN EUR PNG SAS; do
    awk -v pop="$i" '{print $1 "\t" $2 "\t" $3 "\t" $11 "\t" $12 "\t" pop}' \
    $INPUTDIR/all_tag_snps."$i".merged.ALL.0.3_R2_cluster.1KG_phase3_essentials.median_af.bed.extended_LD >> $OUTPUTDIR/ancestral_introgressed_variants_medianextendedLD.bed
done
# sort by chromosome
sort -u $OUTPUTDIR/ancestral_introgressed_variants_medianextendedLD.bed | sort -k 1,1 -k2,2n > $OUTPUTDIR/tmp; mv $OUTPUTDIR/tmp $OUTPUTDIR/ancestral_introgressed_variants_medianextendedLD.bed