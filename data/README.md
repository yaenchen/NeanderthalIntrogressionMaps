## Neanderthal Introgression Maps

This repository contains data and scripts used in the Neanderthal introgression maps comparison. Please see the associated manuscript, "Comparing Neanderthal introgression maps reveals core agreement but substantial heterogeneity."

Unlifted regions from hg19 -> hg38 and T2T are located in each genome build's *unlifted* directory.

For hg19, hg38, and T2T genome builds, we provide:

- all_deserts.bg: contains desert regions across methods, with each column (method) containing either 1 (desert region for that method) or 0 (not a desert regions for that method)

## introgression_maps: Contains introgressed regions and variants. 
NOTE: Files listed directly below are gzipped to save space.

- **all_methods_boolean.bg** contains introgressed regions across methods, with each column (method) containing True or False
- **all_methods.bg** contains introgressed regions across methods, with each column (method) containing an associated score specific to each method, described in the Methods section of the manuscript.
- **84_shared_CEU.bg** contains the introgressed regions for 84 shared 1KG CEU individuals across 6 introgression maps, and is the raw data used to generate Figure 5. 
- **african_introgressed.bg** contains introgressed regions identified by IBDMix, ARGweaver-D, and SARGE. IBDMix regions are based on LWK, GWD, MSL, YRI, ASW, ACB, and ESN from 1KG. ARGweaver-D regions are based on Mandenka_2F and Khomani_San_1F from SGDP. SARGE regions are based on African superpopulation genomes from SGDP (Africa and Africa-MBK, as described in https://www.science.org/doi/10.1126/sciadv.abc0776).


### For each introgression map directory (maps_by_method):
- **METHOD_regions_merged.bed** contains merged *regions*, sorted and merged across populations.
- **METHOD_regions_superpop.bg** contains merged *regions* across superpopulations, with a 1 if it is seen in a population and 0 if not.
- **METHOD_variants.bed** contains *variants* across all populations, separated by column after chrom, start, and end.

For each map with varying amounts of support (maps_by_support):

- **min_num_methods** contains merged regions and variants supported by different numbers of methods. NOTE: Files in this directory only are gzipped to save space.
- **num_methods** contains merged regions and variants (introgressed_regions_n_methods.bed) supported by at least n numbers of methods
