# set up directories
INPUTDIR="/wynton/home/capra/ychen39/introgression_methods/data/sankararaman_2016/summaries/2/neandertal"
OUTPUTDIR="/wynton/home/capra/ychen39/introgression_methods/cleaned/sankararaman_2016"

# add column names
echo -e "Chromosome\tStart\tAlleleScore\tHaplotypeScore\tPopulation" > $OUTPUTDIR/pops.bed
# add data from each population file
for pop in $(cat $INPUTDIR/../../pops); do
    cd $INPUTDIR/"$pop"/summaries/
    FILES=`ls --color=never pred*`
    for file in $FILES; do
        echo $file
        # keep entries if their nean probabilities (column 11) are greater than 0.9
        awk -v pop="$pop" '{ if ($10+0 > 0.5) print $2 "\t" $4 "\t" $10 "\t" $11 "\t" pop}' $file >> $OUTPUTDIR/pops.bed
    done
done

# remove extra file lines except ones starting with numbers (chromosomes), X (for chromosome X), and C (header line)
grep '^[0-9 | X | C]' $OUTPUTDIR/pops.bed > $OUTPUTDIR/neanderthal_introgressed_variants_pops.bed

# set up directories
INPUTDIR="/wynton/home/capra/ychen39/introgression_methods/data/sankararaman_2016/summaries/2/neandertal"
OUTPUTDIR="/wynton/home/capra/ychen39/introgression_methods/cleaned/sankararaman_2016/2"

# add column names
echo -e "Chromosome\tStart\tStop\tHaplotypeScore\tPopulation" > $OUTPUTDIR/fragments.bed

# add data from each population file
for pop in $(cat $INPUTDIR/../../pops); do
    cd $INPUTDIR/"$pop"/summaries/haplotypes
    FILES=`ls --color=never *.haplotypes`
    for file in $FILES; do
        echo $file

        # keep entries if their nean probabilities (column 11) are greater than 0.9
        awk -v pop="$pop" '{print $1 "\t" $3 "\t" $4 "\t" $7 "\t" pop}' $file >> $OUTPUTDIR/fragments.bed
    done
done
grep '^[0-9 | X | C]' $OUTPUTDIR/fragments.bed > $OUTPUTDIR/neanderthal_introgressed_fragments_pops.bed

### Repeat for first directory
# set up directories
INPUTDIR="/wynton/home/capra/ychen39/introgression_methods/data/sankararaman_2016/summaries/1/neandertal"
OUTPUTDIR="/wynton/home/capra/ychen39/introgression_methods/cleaned/sankararaman_2016/1"

mkdir $OUTPUTDIR

# add column names
echo -e "Chromosome\tStart\tAlleleScore\tHaplotypeScore\tPopulation" > $OUTPUTDIR/pops.bed
# add data from each population file
for pop in $(cat $INPUTDIR/../../pops); do
    cd $INPUTDIR/"$pop"/summaries/
    FILES=`ls --color=never pred*`
    for file in $FILES; do
        echo $file
        # keep entries if their nean probabilities (column 11) are greater than 0.9
        awk -v pop="$pop" '{ if ($10+0 > 0.5) print $2 "\t" $4 "\t" $10 "\t" $11 "\t" pop}' $file >> $OUTPUTDIR/pops.bed
    done
done

# remove extra file lines except ones starting with numbers (chromosomes), X (for chromosome X), and C (header line)
grep '^[0-9 | X | C]' $OUTPUTDIR/pops.bed > $OUTPUTDIR/neanderthal_introgressed_variants_pops.bed

# add column names
echo -e "Chromosome\tStart\tStop\tHaplotypeScore\tPopulation" > $OUTPUTDIR/fragments.bed

# add data from each population file
for pop in $(cat $INPUTDIR/../../pops); do
    cd $INPUTDIR/"$pop"/summaries/haplotypes
    FILES=`ls --color=never *.haplotypes`
    for file in $FILES; do
        echo $file
        # keep entries if their nean probabilities (column 11) are greater than 0.9
        awk -v pop="$pop" '{print $1 "\t" $3 "\t" $4 "\t" $7 "\t" pop}' $file >> $OUTPUTDIR/fragments.bed
    done
done
grep '^[0-9 | X | C]' $OUTPUTDIR/fragments.bed > $OUTPUTDIR/neanderthal_introgressed_fragments_pops.bed