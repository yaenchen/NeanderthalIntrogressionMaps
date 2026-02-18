#!/bin/env python

import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import pyBigWig
import pybedtools # if on M1 chip mac while installing bedtools make sure conda config --env --set subdir osx-64 (i was on arm-64 and it could not find bedtools!!!) also install rosetta2 /usr/sbin/softwareupdate --install-rosetta --agree-to-license
from pybedtools import BedTool
import numpy as np
import gzip
import sys
import re
sys.path.append("/wynton/home/capra/ychen39/introgression_methods/introgression_tools/") # MUST ADDS THIS PATH TO IMPORT INTROGRESSION_TOOLS
from introgression_tools import tools
import sys # add scratch dir from bash


# change parent directory to use tools

from pathlib import Path
os.chdir('/wynton/home/capra/ychen39/introgression_methods/introgression_tools')
from introgression_tools import tools

HOMEDIR='/wynton/home/capra/ychen39/'
BEDTOOLSDIR = '/wynton/home/cbi/shared/software/CBI/bedtools2-2.31.1'
# set bedtools path
pybedtools.helpers.set_bedtools_path(path=f'{BEDTOOLSDIR}')

print('scratch dir:', sys.argv[1])


# Read in individual-level files

archaicseeker2_individual=pd.read_csv(f'{HOMEDIR}/introgression_methods/cleaned/introgression_tools/yuan_2021/individual_nean_yuan2021_frag.bed', sep='\t', low_memory=False)
archaicseeker2_individual['ID'] = archaicseeker2_individual['ID'].str.split('_').str[0]
archie_individual=pd.read_csv(f'{HOMEDIR}/introgression_methods/cleaned/introgression_tools/durvasula_2019/individual_nean_durvasula2019_frag.bed', sep='\t', low_memory=False)
dicaladmix_individual=pd.read_csv(f'{HOMEDIR}/introgression_methods/cleaned/introgression_tools/steinruecken_2018/individual_introgressed_steinruecken2018_frag.bed', sep='\t', low_memory=False)
dicaladmix_individual['Chromosome'] = dicaladmix_individual['Chromosome'].str.replace('chr', '')
sankararaman_2014_individual=pd.read_csv(f'{HOMEDIR}/introgression_methods/cleaned/introgression_tools/sankararaman_2014/1kgids_individual_neanderthal_introgressed_fragments_pops.bed', sep='\t', low_memory=False)
vernot_2016_individual=pd.read_csv(f'{HOMEDIR}/introgression_methods/cleaned/introgression_tools/vernot_2016/neanderthal_introgressed_haplotypes_individual.bed', sep='\t', low_memory=False)
admixfrog_individual=pd.read_csv(f'{HOMEDIR}/introgression_methods/cleaned/introgression_tools/iasi_2024/individual_iasi_2024_introgressed_frag.bg', sep='\t', low_memory=False)
ibdmix_2024_individual=pd.read_csv(f'{HOMEDIR}/introgression_methods/cleaned/introgression_tools/li_2024/individual_ibdmix_2024_introgressed_frag.bg', sep='\t', low_memory=False)
hmmix18_individual=pd.read_csv(f'{HOMEDIR}/introgression_methods/cleaned/introgression_tools/skov_2018/individual_nean_skov2018_frag.bed', sep='\t', low_memory=False)

dfs_dict = {'ArchaicSeeker2':archaicseeker2_individual, 
       'ArchIE':archie_individual,
        'DICALADMIX':dicaladmix_individual,
           'Sankararaman14':sankararaman_2014_individual,
           'SStar':vernot_2016_individual,
           'admixfrog':admixfrog_individual,
           'IBDMix24':ibdmix_2024_individual,
           'hmmix18':hmmix18_individual}


# 84 CEU (archie)
heatmap_data=tools.individual_jaccard_comparison(dfs_dict = dfs_dict,
    temp_dir=f'{sys.argv[1]}', 
    output_jaccards_file_path=f'{HOMEDIR}/introgression_methods/cleaned/introgression_tools/individual_comparison/ceu_individual_shuffle.tsv',
    # also save the ibdmix_archaicseeker2_archie_dicaladmix_sstar_sankararaman14_individual file so we can compute lengths distributions across these methods
    # mainly to compare against archie      
    save_shared_individuals_bed=f'{HOMEDIR}/introgression_methods/cleaned/introgression_tools/ceu_individual_shuffle.bg',
    return_heatmap=True,
    autosomes_only=True,
    x_only=False,
    genome_file = f'{HOMEDIR}/introgression_methods/data/human_hg19_gaps_removed_common_chr_autosomes_no_prefix.bed',)

# # 1KG (no archie) - should be 271 CEU and CHBS individuals
# dfs_dict = {'IBDMix20':ibdmix_individual, 
#        'ArchaicSeeker2':archaicseeker2_individual, 
#         'DICALADMIX':dicaladmix_individual,
#            'Sankararaman14':sankararaman_2014_individual,
#            'SStar':vernot_2016_individual,
#            'admixfrog':admixfrog_individual,
#            'IBDMix24':ibdmix_2024_individual,
#            'hmmix18':hmmix18_individual}

# heatmap_data=tools.individual_jaccard_comparison(dfs_dict = dfs_dict,
#     temp_dir=f'{sys.argv[1]}', 
#     output_jaccards_file_path=f'{HOMEDIR}/introgression_methods/cleaned/introgression_tools/individual_comparison/271CEUCHBS.tsv',
#     # also save the ibdmix_archaicseeker2_archie_dicaladmix_sstar_sankararaman14_individual file so we can compute lengths distributions across these methods
#     # mainly to compare against archie      
#     save_shared_individuals_bed=f'{HOMEDIR}/introgression_methods/cleaned/introgression_tools/individual_comparison/271CEUCHBS.bg',
#     return_heatmap=True,
#     autosomes_only=True,
#     x_only=False,
#     paired_shuffled_match=True,
#     genome_file = f'{HOMEDIR}/introgression_methods/data/human_hg19_gaps_removed_common_chr_autosomes_no_prefix.bed',)

# 1 SGDP Papuan individual
