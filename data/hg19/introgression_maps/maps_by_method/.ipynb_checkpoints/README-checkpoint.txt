## Uploaded by Yaen Chen, 12/2/2024

## Contains introgressed regions and variants across methods

For each introgression map:
    - sorted_union_merged_chr_prefix.bg contains merged *regions*, sorted and merged across populations with 'chr' prefix in the first column -- Can be used with make_annots for ldsc
    - union_merged.bg contains merged *regions*, with a 1 if it is seen in a population and 0 if not.
    - neanderthal_introgressed_variants_chr_prefix.bg contains *variants* across all populations, merged and with 'chr' prefix in the first column -- Can be used with make_annots for ldsc