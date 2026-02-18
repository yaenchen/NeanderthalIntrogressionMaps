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
import glob
import pickle
import time


# change parent directory to use tools
from pathlib import Path
os.chdir('/wynton/home/capra/ychen39/introgression_methods/introgression_tools')
from introgression_tools import tools
HOMEDIR='/wynton/home/capra/ychen39/'
BEDTOOLSDIR = '/wynton/home/cbi/shared/software/CBI/bedtools2-2.31.1/bin'
# set bedtools path
pybedtools.helpers.set_bedtools_path(path=f'{BEDTOOLSDIR}')

# autosomes only
autosomes = 'True'

# set scratch dir
print('scratch dir:', sys.argv[1])
scratch_dir = sys.argv[1]

# number of shuffles
num_shuffles = int(sys.argv[2])

# make output dir
output_dir = sys.argv[3]
os.system(f"mkdir {output_dir}")

# directories for each method's data
methods_dirs = ['yuan_2021', 'durvasula_2019', 'steinruecken_2018', 'skov_2020', 'schaefer_2021', 'hubisz_2020', 'browning_2018_updated_filtering', 
                'vernot_2016', 'sankararaman_2014', 'sankararaman_2016_1', 'sankararaman_2016_2', 'li_2024', 'skov_2018', 'iasi_2024']

# renaming dict
rename_dict = {'ibdmix20': 'IBDMix20', 'argweaverd': 'ARGweaver-D', 'archaicseeker': 'ArchaicSeeker2', 'archie': 'ArchIE',
                                   'sprime': 'Sprime', 'sarge': 'SARGE', 'vernot_2016':'S*', 'sankararaman_2014': 'CRF14',
                                 'sankararaman_2016_1':'CRF16 (1)', 'sankararaman_2016_2':'CRF16 (2)', 'skov_2020':'Skov20',
                                 'steinruecken_2018':'DICAL-ADMIX', 'vernot_2016_extendedLD': 'S* (extended LD)', 
                                 'vernot_2016_medianextendedLD': 'S* (median extended LD)', 'ibdmix24':'IBDMix', 'iasi_2024':'admixfrog', 'skov_2018':'hmmix18'} # added after reviews


# read in all the data
# read in all the data
# make a dictionary to get from Illumina_ID to SGDP_ID
SGDP_metadata=pd.read_csv(f'{HOMEDIR}/introgression_methods/data/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt', sep='\t', low_memory=False, encoding='unicode_escape')
SGDP_metadata['SGDP_ID'] = SGDP_metadata['SGDP_ID'].apply(lambda x: '_'.join(x.split('_')[1:]))
SGDP_metadata_dict = dict(zip(SGDP_metadata.Illumina_ID, SGDP_metadata.SGDP_ID))

SARGE_individual=pd.read_csv(f'{HOMEDIR}/introgression_methods/data/sarge_schaefer2021/sarge_admixed_haps_hs37d5.bed', sep='\t', low_memory=False)
SARGE_individual.rename(columns={'#chrom': 'Chromosome',
                                'start': 'Start',
                                'end': 'End',
                                'len': 'length',
                                'score': 'Score'}, inplace=True)

# remove prefix and suffix
SARGE_individual['ID'] = SARGE_individual['hap'].apply(lambda x: '_'.join(x.split('_')[1:]))
SARGE_individual['ID'] = SARGE_individual['ID'].apply(lambda x: '-'.join(x.split('-')[:2]))

# make a dictionary to get from Illumina_ID to SGDP_ID
SGDP_metadata=pd.read_csv(f'{HOMEDIR}/introgression_methods/data/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt', sep='\t', low_memory=False, encoding='unicode_escape')
SGDP_metadata['SGDP_ID'] = SGDP_metadata['SGDP_ID'].apply(lambda x: '_'.join(x.split('_')[1:]))
SGDP_metadata_dict = dict(zip(SGDP_metadata.Illumina_ID, SGDP_metadata.SGDP_ID))

admixfrog_individual=pd.read_csv(f'{HOMEDIR}/introgression_methods/cleaned/introgression_tools/iasi_2024/individual_iasi_2024_introgressed_frag.bg', sep='\t', low_memory=False)
# remove prefix and suffix
admixfrog_individual['ID'] = admixfrog_individual['ID'].apply(lambda x: '_'.join(x.split('_')[1:]))
admixfrog_individual['ID'] = admixfrog_individual['ID'].apply(lambda x: '-'.join(x.split('-')[:2]))

sankararaman16_1_individual=pd.read_csv(f"{HOMEDIR}/introgression_methods/cleaned/introgression_tools/sankararaman_2016_1/sgdp_individual_neanderthal_introgressed_fragments_pops.bed", sep='\t', low_memory=False)
sankararaman16_1_individual['ID'] = sankararaman16_1_individual['ID'].apply(lambda x: '_'.join(x.split('_')[1:]))
sankararaman16_1_individual['ID'] = sankararaman16_1_individual['ID'].apply(lambda x: '-'.join(x.split('-')[:2]))
sankararaman16_2_individual=pd.read_csv(f"{HOMEDIR}/introgression_methods/cleaned/introgression_tools/sankararaman_2016_2/sgdp_individual_neanderthal_introgressed_fragments_pops.bed", sep='\t', low_memory=False)
sankararaman16_2_individual['ID'] = sankararaman16_2_individual['ID'].apply(lambda x: '_'.join(x.split('_')[1:]))
sankararaman16_2_individual['ID'] = sankararaman16_2_individual['ID'].apply(lambda x: '-'.join(x.split('-')[:2]))

dfs_dict = {'SARGE':SARGE_individual, 
            'CRF16_1':sankararaman16_1_individual,
            'CRF16_2':sankararaman16_2_individual,
            'admixfrog':admixfrog_individual,
            }

del SARGE_individual, sankararaman16_1_individual, sankararaman16_2_individual, admixfrog_individual

num_algorithms = len(dfs_dict)

algorithms = list(dfs_dict.keys())
# group by ID and count unique values in 'Value'
individual_concat = pd.concat([df.assign(Algorithm=k) for k,df in dfs_dict.items()])
id_value_counts = individual_concat.groupby('ID')['Algorithm'].nunique()
shared_individuals = id_value_counts[id_value_counts == num_algorithms].index.tolist()    

for person in shared_individuals:
    for algorithm_name in algorithms:
        # create a bedgraph for that individual
        if autosomes:
            person_temp = individual_concat[(individual_concat['ID']==person) & (individual_concat['Algorithm']==algorithm_name) & (individual_concat['Chromosome'] != 'X')][['Chromosome', 'Start', 'End', 'Score']]
            person_temp = person_temp[person_temp['Chromosome']!='X']
        else:
            person_temp = individual_concat[(individual_concat['ID']==person) & (individual_concat['Algorithm']==algorithm_name) & (individual_concat['Chromosome'] == 'X')][['Chromosome', 'Start', 'End', 'Score']]

        person_temp = individual_concat[(individual_concat['ID']==person) & (individual_concat['Algorithm']==algorithm_name)][['Chromosome', 'Start', 'End', 'Score']]
        # temporarily save the dfs to run bedtools on them
        person_temp.to_csv(f'{scratch_dir}/{algorithm_name}_{person}', sep='\t', index=False, header=False)
        # sort
        os.system(f"sort -k1,1 -k2,2n {scratch_dir}/{algorithm_name}_{person} -o {scratch_dir}/{algorithm_name}_{person}")
        # remove X if autosomes using sed inplace - just in case it was missed earlier
        if autosomes:
            os.system(f"sed -i '/^X/d' {scratch_dir}/{algorithm_name}_{person}")

        for i in range(int(num_shuffles)):
            # create shuffled bedgraph
            os.system(f'{BEDTOOLSDIR}/bedtools shuffle -i {scratch_dir}/{algorithm_name}_{person} -g {HOMEDIR}/introgression_methods/data/human_hg19_gaps_removed_common_chr_autosomes_no_prefix.bed > {scratch_dir}/{algorithm_name}_Shuffled_{person}_{i}.bg')
            # sort
            os.system(f"sort -k1,1 -k2,2n {scratch_dir}/{algorithm_name}_Shuffled_{person}_{i}.bg -o {scratch_dir}/{algorithm_name}_Shuffled_{person}_{i}_sort.bg")
            # merge
            os.system(f'{BEDTOOLSDIR}/bedtools merge -i {scratch_dir}/{algorithm_name}_Shuffled_{person}_{i}_sort.bg -c 4 -o max > {scratch_dir}/{algorithm_name}_Shuffled_{person}_{i}.bg')

#read in the shuffles and compute jaccard similarities
#list to store final heatmap data
final_jaccards_shuffled = []
final_jaccards_observed = []
final_z_scores = []
final_empirical_p_values = []
person_ids = [] # this one also stores IDs, in the order that they are in the final tsv file
# for each individual
for person in shared_individuals:
    # progress print
    print('Processing individual:', person)
    # save jaccard similarities, z scores, p values for each person
    jaccard_similarities={}
    z_scores = pd.DataFrame(index=algorithms, columns=algorithms)
    p_values = pd.DataFrame(index=algorithms, columns=algorithms)
    # for each shuffle
    for i in range(num_shuffles):
        # get input files for each shuffle round
        input_files = glob.glob(f'{scratch_dir}/*_Shuffled_{person}_{i}.bg')
        # separate by space
        input_files = ' '.join([file for file in input_files])
        # to get names, keep only the method names from input_files
        if scratch_dir.endswith('/'):
            algorithms_str = input_files.replace(f'{scratch_dir}', '').replace(f'_Shuffled_{person}_{i}.bg', '')
        else:
            algorithms_str = input_files.replace(f'{scratch_dir}/', '').replace(f'_Shuffled_{person}_{i}.bg', '')
        overlap = tools.generate_overlap(
                input_files_string=input_files,
                input_files_names_str=f'{algorithms_str}',
                                        output_file_path=f'{scratch_dir}/shuffle_{person}_{i}',
                                        boolean=True)

        # compute jaccard similarities of shuffles
        for cat1 in algorithms:
            for cat2 in algorithms:
                jaccard_similarities[(cat1, cat2)] = tools.jaccard_similarity(overlap, cat1, cat2, 'length')

        # create heatmap df
        heatmap_data_shuffled = pd.DataFrame(index=algorithms, columns=algorithms)
        # save jaccards to heatmap df
        for (cat1, cat2), similarity in jaccard_similarities.items():
            heatmap_data_shuffled.loc[cat1, cat2] = similarity

        # save jaccards from each shuffle
        final_jaccards_shuffled.append(heatmap_data_shuffled)


    # NOW FOR THE OBSERVED DATA
    observed_jaccard_similarities = {}
    # compute raw jaccards
    input_files = glob.glob(f'{scratch_dir}/*_{person}')
    # remove observed
    input_files = [file for file in input_files if 'observed' not in file]
    # separate by space
    input_files = ' '.join([file for file in input_files])
    # to get names, keep only the method names from input_files
    # remove observed
    if scratch_dir.endswith('/'):
        algorithms_str = input_files.replace(f'{scratch_dir}', '').replace(f'_{person}', '')
    else:
        algorithms_str = input_files.replace(f'{scratch_dir}/', '').replace(f'_{person}', '')
    overlap = tools.generate_overlap(
            input_files_string=input_files,
            input_files_names_str=f'{algorithms_str}',
                                    output_file_path=f'{scratch_dir}/observed_{person}',
                                    boolean=True)   
    # compute jaccard similarities of observed
    for cat1 in algorithms:
        for cat2 in algorithms:
            observed_jaccard_similarities[(cat1, cat2)] = tools.jaccard_similarity(overlap, cat1, cat2, 'length')

    # create heatmap df for observed
    heatmap_data_observed = pd.DataFrame(index=algorithms, columns=algorithms)
    # save jaccards to heatmap df
    for (cat1, cat2), similarity in observed_jaccard_similarities.items():
        heatmap_data_observed.loc[cat1, cat2] = similarity

    #calculate z scores, empirical p values for each observed value against shuffled distribution
    for cat1 in algorithms:
        for cat2 in algorithms:
            observed_value = heatmap_data_observed.loc[cat1, cat2]
            # get shuffled values
            shuffled_values = [df.loc[cat1, cat2] for df in final_jaccards_shuffled]
            # calculate z score
            if len(shuffled_values) > 0:
                mean = np.mean(shuffled_values)
                std = np.std(shuffled_values)
                if std != 0:
                    z_scores.loc[cat1, cat2] = (observed_value - mean) / std
                else:
                    z_scores.loc[cat1, cat2] = 0
            else:
                z_scores.loc[cat1, cat2] = 0
            # calculate empirical p value
            count = sum(1 for value in shuffled_values if value >= observed_value)
            p_values.loc[cat1, cat2] = (count + 1) / (len(shuffled_values) + 1)

    # save z scores and p values across people
    final_z_scores.append(z_scores)
    final_empirical_p_values.append(p_values)
    final_jaccards_observed.append(heatmap_data_observed)
    # add ID to running list to add to final dataframe
    person_ids.append(person)

# after computing z scores and p values for each individual, save to a final dataframe with renamed columns
print('Saving final dataframes...')
for df_name, df in {'Jaccards': final_jaccards_observed, 'Empirical_P_Values': final_empirical_p_values, 'Z_Scores': final_z_scores}.items():
    # flatten every heatmap dataframe -> single column
    flattened_columns = []
    for i, person_df in enumerate(df):
        # melt dataframe to long format
        melted_df = person_df.reset_index().melt(id_vars='index', var_name='column', value_name=f'{df_name}')
        # create combined identifier for the pairs (method1:method2)
        melted_df['Pair'] = melted_df['index'] + ':' + melted_df['column']
        # set 'pair' column as index
        melted_df.set_index('Pair', inplace=True)
        # add corresponding person ID
        melted_df['ID'] = person_ids[i]  
        # append melted DataFrame
        flattened_columns.append(melted_df[['ID', f'{df_name}']])

        # combine all into final dataframe
        flat_df = pd.concat(flattened_columns).reset_index().rename(columns={'index': 'Pair'})
        # save final dataframe
        flat_df.to_csv(f'{output_dir}/{df_name}', sep='\t', index=False)

# delete temp files
# os.system(f'rm {scratch_dir}/*')

# save heatmap_data_observed to pickle file
with open(f'{output_dir}/heatmap_data_observed.pkl', 'wb') as f:
    pickle.dump(final_jaccards_observed, f)