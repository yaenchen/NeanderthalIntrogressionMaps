#!/bin/env python

import gzip
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import seaborn as sns
import pyranges as pr
from fuc import pybed
import pyBigWig
import os
import pybedtools
BEDTOOLSDIR = '/wynton/home/cbi/shared/software/CBI/bedtools2-2.31.1/bin'
pybedtools.helpers.set_bedtools_path(path=f'{BEDTOOLSDIR}')
import tarfile
import io
import re
import sys

HOMEDIR = '/wynton/home/capra/ychen39/'
DIR = '/wynton/group/capra/data/introgression_maps/steinruecken_2018/CalledTracts_March2018/'
PROB_DIR = f'{HOMEDIR}/introgression_methods/data/steinruecken_2018/posterior_probabilities'

def collapse_columns(df,col):
    cols = df.columns.values.tolist()
    cols.remove(col)
    df = df.groupby(cols)[col].apply(';'.join).reset_index()
    return df


# Merge all the files and save as a single bed file
groups = ['CEU', 'CHBS']

steinruecken = pd.DataFrame()

for i in groups:
    # get the tar.gz file passed in
    j = str(sys.argv[1:][0])
    print(j)
    # loop through each population at a time
    if j.startswith(i):
        f = DIR + j
        prob_dir = PROB_DIR + '/' + j.replace('_beds', '')
        tar = tarfile.open(f, 'r:gz')
        bed_files = [f.name for f in tar.getmembers() if f.name.endswith('.bed')]
        
        for bed_file in bed_files:
            # get the person ID from the name of the bed file
            id = bed_file.split(".")[1].replace('/', '')
            # extract bed contents
            bed_contents = tar.extractfile(bed_file).read()
            # write bed contents to csv
            df_temp = pd.read_csv(io.BytesIO(bed_contents), sep='\t', header=None)
            df_temp['Population'] = i
            df_temp['ID'] = id
            # get the corresponding probabilities
            # get the file name without ".bed" suffix
            prob_tar_file = PROB_DIR + '/' + bed_file.split("/")[0] + '.tar.gz'
            prob_file = bed_file.split("/")[0] +'/./' + bed_file.split("/")[2].replace('.bed', '.filtered')
            pos_file = bed_file.split("/")[0] +'/./' + bed_file.split("/")[0].split("_")[2] + '.pos'
            prob_tar = tarfile.open(prob_tar_file, 'r:gz')
            prob_contents = prob_tar.extractfile(prob_file).read()
            pos_contents = prob_tar.extractfile(pos_file).read()
            prob_df_temp = pd.read_csv(io.BytesIO(prob_contents), sep='\t', header=None)
            pos_df_temp = pd.read_csv(io.BytesIO(pos_contents), sep='\t', header=None)
            # get the positions and probailities
            position_probs = pd.concat([pos_df_temp.reset_index(drop=True), prob_df_temp.reset_index(drop=True)], axis=1)
            position_probs.columns=['End','Score']
            df_temp.columns=['Chromosome','Start', 'End', 'Ancestry', 'ID']
            df_temp[['Chromosome', 'Ancestry', 'ID']] = df_temp[['Chromosome', 'Ancestry', 'ID']].astype(str)
            # merge the data only if the positions appear in both
            final = df_temp.merge(position_probs, on='End')
            steinruecken = pd.concat([steinruecken, final])

OUTPUT_FILE = f'{HOMEDIR}/introgression_methods/cleaned/steinruecken_2018/round2/{j}_frag.bed'
steinruecken.rename(columns={0:'Chr',1:'Start',2:'End', 3:'Ancestry', 4:'ID', 5:'Score'}, inplace=True)
steinruecken.to_csv(OUTPUT_FILE, sep='\t', index=False)