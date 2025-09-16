### Generates variant and haplotype level files

# set up directories
# local
INPUTDIR="/Users/yaenchen/Library/CloudStorage/Box-Box/wynton_backup/introgression_methods/data/sankararaman_2014_copy/sankararaman_2014"
OUTPUTDIR="/Users/yaenchen/Library/CloudStorage/Box-Box/wynton_backup/introgression_methods/cleaned/sankararaman_2014_copy/"
# wynton
#INPUTDIR="/wynton/home/capra/ychen39/introgression_methods/data/sankararaman_2014_copy/sankararaman_2014"
#OUTPUTDIR="/wynton/home/capra/ychen39/introgression_methods/cleaned/sankararaman_2014/"

mkdir $OUTPUTDIR

# for each population files (not extended LD)
# add column names
echo -e "Chromosome\tStart\tAlleleScore\tHaplotypeScore\tAncestry" > $OUTPUTDIR/neanderthal_introgressed_variants_pops.bed
# add data from each population file
for pop in $(cat $INPUTDIR/pops.EUR-ASN); do
    cd $INPUTDIR/"$pop".hapmap/summaries/
    FILES=`ls --color=never *.gz`
    for file in $FILES; do
        echo $file
        # keep entries if their nean probabilities (column 11) are greater than 0.9
        gunzip -c $file | awk -v pop="$pop" '{ if ($10+0 > 0.9) print $2 "\t" $4 "\t" $10 "\t" $11 "\t" pop}' >> $OUTPUTDIR/pops.bed
    done
done
# remove extra file lines except ones starting with numbers (chromosomes), X (for chromosome X), and C (header line)
grep '^[0-9 | X | C]' $OUTPUTDIR/pops.bed > $OUTPUTDIR/neanderthal_introgressed_variants_pops.bed

# set up directories
# local
INPUTDIR="/Users/yaenchen/Library/CloudStorage/Box-Box/wynton_backup/introgression_methods/data/sankararaman_2014_copy/sankararaman_2014"
OUTPUTDIR="/Users/yaenchen/Library/CloudStorage/Box-Box/wynton_backup/introgression_methods/cleaned/sankararaman_2014_copy/"
# wynton
#INPUTDIR="/wynton/home/capra/ychen39/introgression_methods_updated/data/sankararaman_2014_copy/sankararaman_2014"
#OUTPUTDIR="/wynton/home/capra/ychen39/introgression_methods_updated/cleaned/sankararaman_2014/"

# for each population files (not extended LD)
# add column names
echo -e "Chromosome\tStart\tEnd\tScore\tAncestry" > $OUTPUTDIR/neanderthal_introgressed_fragments_pops.bed
# add data from each population file
for pop in $(cat $INPUTDIR/pops.EUR-ASN); do
    cd $INPUTDIR/"$pop".hapmap/summaries/haplotypes
    FILES=`ls --color=never *.gz`
    for file in $FILES; do
        echo $file
        # keep entries if their nean probabilities (column 11) are greater than 0.9
        gunzip -c $file | awk -v pop="$pop" '{print $1 "\t" $3 "\t" $4 "\t" $7 "\t" pop}' >> $OUTPUTDIR/fragments.bed
    done
done
 # remove extra file lines except ones starting with numbers (chromosomes), X (for chromosome X), and C (header line)
grep '^[0-9 | X | C]' $OUTPUTDIR/fragments.bed > $OUTPUTDIR/neanderthal_introgressed_fragments_pops.bed

