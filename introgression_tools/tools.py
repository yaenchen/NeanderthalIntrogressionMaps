# Imports
import gzip
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pybedtools import BedTool
import os
import pybedtools # if on M1 chip mac while installing bedtools make sure conda config --env --set subdir osx-64 (i was on arm-64 and it could not find bedtools!!!) also install rosetta2 /usr/sbin/softwareupdate --install-rosetta --agree-to-license
import tarfile
import io
import sys
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
bedops_dir='/wynton/home/cbi/shared/software/CBI/bedops-2.4.41/bin/'

# manually changed this path
def get_bedtools_path_from_user():
    while True:
        #bedtools_path = "/opt/local/bin/"
        bedtools_path = '/wynton/home/cbi/shared/software/CBI/bedtools2-2.31.1/bin/'
        #bedtools_path = input("Enter the bedtools directory path: ")
        if os.path.exists(bedtools_path):
            return bedtools_path
        else:
            print("Invalid directory path. Please try again.")

# set bedtools path for pybedtools
BEDTOOLSDIR = get_bedtools_path_from_user()
pybedtools.helpers.set_bedtools_path(path=f'{BEDTOOLSDIR}')

# initialize population codes -> superpopulation dictionary (1000 Genomes)
ancestry_dict = {'EAS': ['EAS', 'CDX', 'CHB', 'CHBS', 'CHS', 'JPT', 'KHV', 'eastasia', 'EastAsia'],
                'SAS': ['SAS', 'BEB', 'GIH', 'ITU', 'PJL', 'STU', 'southasia', 'SouthAsia'],
                'EUR': ['EUR', 'CEU', 'FIN', 'GBR', 'IBS', 'TSI', 'westeurasia', 'WestEurasia'],
                'AMR': ['AMR', 'PEL', 'CLM', 'MXL', 'PUR', 'america', 'America'],
                'Oceania': ['Papuans', 'Papuan', 'oceania', 'Oceania', 'PNG'],
                'CentralAsia': ['CentralAsia', 'centralasia', 'CentralAsiaSiberia'],
                'AFR': ['AFR', 'Africa', 'Africa2']}

# dictionary to remap algorithm names
rename_dict = {'ibdmix': 'IBDMix', 'argweaverd': 'ARGweaver-D', 'archaicseeker': 'ArchaicSeeker2', 'archie': 'ArchIE',
                                   'sprime': 'Sprime', 'sarge': 'SARGE', 'vernot_2016':'S*', 'sankararaman_2014': 'Sankararaman14',
                                 'sankararaman_2016_1':'Sankararaman16 (1)', 'sankararaman_2016_2':'Sankararaman16 (2)', 'skov_2020':'Skov20',
                                 'steinruecken_2018':'DICAL-ADMIX'}

# function to map population to superpopulation
def map_status(anc):
    for key, values in ancestry_dict.items():
        if anc in values:
            return key

def save_bed_only(df=False, full_file_path=False, output_file_path=False, colnames=False, prefix=False, save_with_header=False, stratify=False, output_dir=False, ancestry=False):
    '''
    take an input file with multiple columns and only keep the bed contents (chromosome, start, end columns)

    Args:
        df: pandas dataframe to input if not a file name
        full_file_path: full path of the input file
        output_file_path: full path of the output file to save as a bed only file
        colnames: list of column names, or False if already in the data
        prefix: True if you want to keep the 'chr' prefix, False if not
        save_with_header: True if you want to save the header, otherwise default no header saved
        stratify: str for a column to stratify by. Lets you stratify files by categories (such as population), and writes each bed file to a separate file based on category + suffix). must use output_dir
        ancestry: if stratifying by ancestry, use ancestry=True to remap population codes to ancestral group codes (e.g. CEU -> EUR)

    NOTE: sorted using bedops, bedtools was giving me trouble? also it lets me use a unique parameter to report only unique regions
    '''
    # check if it's a dataframe, if not, read as a file
    if isinstance(df, pd.DataFrame) != True and full_file_path != False:
        if colnames:
            df=pd.read_csv(full_file_path, names=colnames, sep='\t', low_memory=False, index_col=False)
        else:
            df=pd.read_csv(full_file_path, sep='\t', low_memory=False, index_col=False)
        
        
    # rename columns based on different possible naming
    df.rename(columns=lambda col: 'chrom' if col in ['chromosome', 'Chromosome', 'Chrom', 'Chr', 'chr'] else col, inplace=True)
    df.rename(columns=lambda col: 'start' if col in ['Start'] else col, inplace=True)
    df.rename(columns=lambda col: 'end' if col in ['End'] else col, inplace=True)
    
    # if missing end position for snps
    if 'end' not in df.columns:
        df['end'] = df['start'].astype(int) + 1
    elif 'start' not in df.columns:
        df['start'] = df['end'] - 1

    df['chrom'] = df['chrom'].astype(str)  # Convert to string for consistent manipulation

    # handle prefix based on the user input
    if prefix==True:
        if df['chrom'].str.contains('chr').all()==False:
            df['chrom'] = 'chr' + df['chrom']
    elif prefix==False:
        df['chrom'] = df['chrom'].str.replace('chr', '')

    # remove duplicates
    df=df.drop_duplicates()
    
    # if a stratification category is specified:
    if stratify:
        # make output dir
        os.system(f'mkdir {output_dir}')
        # subset df by cat
        df=df[['chrom', 'start', 'end', stratify]]
        
        if ancestry == True:
            df[stratify]=df[stratify].apply(map_status)
        
        for category in df[stratify].unique():
            subset=df[df[stratify]==category]
            
            subset=subset[['chrom', 'start', 'end']]
            
            subset=subset.drop_duplicates()
            
            subset.to_csv(f'{output_dir}/{category}_unsorted', header=save_with_header, sep='\t', index=False)
            
            if prefix==True:
                # sort file and keep unique snps
                os.system(f'{bedops_dir}/sort-bed --unique {output_dir}/{category}_unsorted > {output_dir}/{category}_chr_prefix.bed')
                os.system(f'rm {output_dir}/{category}_unsorted')
                      
            else:
                # sort file and keep unique snps
                os.system(f'{bedops_dir}/sort-bed --unique {output_dir}/{category}_unsorted > {output_dir}/{category}.bed')
                os.system(f'rm {output_dir}/{category}_unsorted')
                      
    # if no stratification:        
    else:    
        df=df[['chrom', 'start', 'end']]
        
        df.to_csv(f'{output_file_path}_unsorted', header=save_with_header, sep='\t', index=False)
        
        # sort file and keep unique snps
        os.system(f'{bedops_dir}/sort-bed --unique {output_file_path}_unsorted > {output_file_path}_sorted')
        os.system(f'rm {output_file_path}_unsorted')

        # merge in case they are regions
        os.system(f'{BEDTOOLSDIR}/bedtools merge -i {output_file_path}_sorted > {output_file_path}_merged')
        os.system(f'rm {output_file_path}_sorted')

        # sort file and keep unique snps
        os.system(f'sort -k1,1 -k2,2n {output_file_path}_merged | uniq > {output_file_path}')
        os.system(f'rm {output_file_path}_merged')        
        
def save_bedgraph(df, save_path):
    """
    This method saves a .bed file into bedgraph format
    Args:
        df: pandas df to save
        save_path: full file path to save output
    """

    for ancestry in df['Ancestry'].unique():
        if ancestry != (None or 'None'):
            print(ancestry)
            bedgraph = df[(df['Ancestry'] == ancestry)]
            bedgraph = bedgraph[['Chromosome', 'Start', 'End', 'Score']]
            if bedgraph['Chromosome'] is str:
                if bedgraph['Chromosome'].str.contains('chr').any():
                    bedgraph['Chromosome'] = bedgraph['Chromosome'].str.replace('chr', '')
            bedgraph.to_csv(f'{save_path}/{ancestry}.bg', sep ='\t', index=False, header=False)
            # sort bedgraph files by position
            os.system(f'{BEDTOOLSDIR}/bedtools sort -i {save_path}/{ancestry}.bg > {save_path}/sorted_{ancestry}.bg')
            # merge overlapping entries, keeping the highest score
            os.system(f'{BEDTOOLSDIR}/bedtools merge -i {save_path}/sorted_{ancestry}.bg -c 4 -o max > {save_path}/merged_{ancestry}.bg')
            print('merged')

def save_to_bed(df, save_path):
    """
    This method saves a pandas DataFrame and converts it into a .bed file

    df: pandas dataframe to be saved
    save_path: string, where you want the .bed file to be saved
    """
    # get headers
    headers = list(df.columns.values)
    headers = [item.strip() for item in headers]
    header_string = '\t'.join(headers)
    #read into PyRanges format
    df_bed = pr.PyRanges(df)
    df_bed = df_bed.df
    # save to bed file
    df_bed.to_file(save_path)
    # add a header to the file
    os.system(f'echo "{header_string}" | cat - "{save_path}" > temp && mv temp "{save_path}"')

def chen_2020_raw_to_bed(input_path, output_dir, output_file_name, population=False, genomes_source=False, individual=False):
    """
    This method reads the raw output from chen_2020 (ibdmix) and cleans it to .bed format with overlapping regions collapsed into one entry, either at the individual or population-levels.
    
    Args:
        input_path: path to .txt of the raw ibdmix file
        output_dir: output directory to place the cleaned file 
        output_file_name: output file name to place in the output_dir
        genomes_source: the project that the genomes are from (usually 1000 Genomes Project (1KG) or Simons Genome Diversity Project (SGDP))
        population: true if to be cleaned at the population level, input_path file should already be cleaned at the individual level
        individual: true if only to be cleaned at the individual level
    """
    df = pd.read_csv(input_path, sep='\t')
    # make output directory
    os.system(f'mkdir {output_dir}')
    # make directory to store bedgraphs
    os.system(f'mkdir {output_dir}/bedgraphs')
    # get final save path
    save_path = f'{output_dir}/{output_file_name}'
    if individual == True:
        df = pd.read_csv(input_path, sep='\t', skiprows=1, names=['Chromosome', 'Start', 'End', 'LOD', 'maxLOD', 'size', 'Population', 'Ancestry', 'ID'])
        # add score_name and fill with value "LOD"
        df['ScoreName'] = 'LOD'
        # add IntrogressionSource (all Neanderthal)
        df['IntrogressionSource'] = 'Neanderthal'
        # add GenomesSource (all 1KG)
        df['Method'] = 'IBDmix, Chen 2020'
        df['GenomesSource'] = '1KG Phase 3'
        # rename LOD to score and change formatting of bed file so we can use pybed to save to .bed
        df.rename(columns={'LOD':'Score'}, inplace=True)
        # drop maxLOD and size columns
        df.drop(['maxLOD', 'size'], axis=1, inplace=True)
        # save to bed file
        save_to_bed(df, f'{save_path}')
    if population == True:
        # convert data to bedgraph format
        save_bedgraph(df, output_dir+'/bedgraphs')
        # get ancestries as a list and string separated by a space
        ancestries_list = list(df['Ancestry'].unique())
        # get list of bedgraph files to be merged in the next step, format is full path separated by " \"
        merged_bg_files = " \\".join([f'{output_dir}'+'/bedgraphs/'+file for file in os.listdir(f'{output_dir}/bedgraphs') if file.startswith("merged_")])
        # get list of ancestries to be used in unionbedg (must be in the same order as the input files)
        ancestries_str = re.sub(f'{output_dir}'+'/bedgraphs/'+r'(merged_|\.bg)', '', merged_bg_files)
        ancestries_str = re.sub(r'(|\.bg)', '', ancestries_str)
        ancestries_str = re.sub(re.escape('\\'), "", ancestries_str)
        # create union file at the population level
        os.system(f'{BEDTOOLSDIR}/bedtools unionbedg -i {merged_bg_files} -header -names {ancestries_str} > {output_dir}/bedgraphs/union_merged.bg')
        # delete intermediate files, keeping only merged files at the ancestral level
        os.system(f'cd {output_dir}/bedgraphs; rm -f sorted*.bg')
        # read in generated union file
        ibdmix_final = pd.read_csv(f'{output_dir}/bedgraphs/union_merged.bg', sep='\t')

        # keep only the highest reported score, collapse into one column
        ibdmix_final['highest_score'] = ibdmix_final[ancestries_list].max(axis=1)
        # drop ancestries
        ibdmix_final = ibdmix_final.drop(ancestries_list, axis=1)
        # save as a bedgraph
        ibdmix_final.to_csv(f'{output_dir}/bedgraphs/union_merged_highest_score.bg', sep='\t', index=False)

        # add metadata
        # add score_name and fill with value "LOD"
        ibdmix_final['ScoreName'] = 'LOD'
        # add IntrogressionSource (all Neanderthal, IBDMix only detects neanderthal introgression)
        ibdmix_final['IntrogressionSource'] = 'Neanderthal'
        # add method name
        ibdmix_final['Method'] = 'IBDmix, Chen 2020'
        # add genomes source (user-inputted parameter)
        if genomes_source:
            ibdmix_final['GenomesSource'] = genomes_source
        else: # assume original input from ibdmix paper, which uses 1KG Phase 3 genomes
            ibdmix_final['GenomesSource'] = '1KG Phase 3'
        # save final file
        ibdmix_final.to_csv(f'{output_dir}/{output_file_name}', sep='\t', index=False)
        
        # completed
        print(f'Output .bed file is located at {save_path}', sep='')

def hubisz_2020_raw_to_bed(output_dir, output_file_name, input_path=None, ancestry=None, ID=None, archaic=None, individual_files=None, genomes_source='SGDP', population=False, individual=False):
    """
    If individual=True, This function takes in an input .bed file (such as ooaM1A/Mandenka_2F.bed) and cleans it into a standardized output
    if population=True and individual_files are provided, it will concatenate the files and generate a population-level file
     
    Args:
        individual_files: list of individual-level file paths
        output_dir: output directory to place the cleaned file, will contain cleaned population/chromosome-level files if individual=True
        output_file_name: output file name to place in the output_dir
        input_dir_list: list of input directories (to be used with individual=True)
        prob_dir: posterior probabilities directory
        genomes_source: the project that the genomes are from (usually 1000 Genomes Project (1KG) or Simons Genome Diversity Project (SGDP))
        file_path = input file containing individual=level introgressed regions (to be used with population=True)
        population: true if to be cleaned at the population level, input_path file should already be cleaned at the individual level
        individual: true if only to be cleaned at the individual level, then must include:
            input_path, output_dir, output_file_name, ancestry, ID, archaic
            returns standardized individual-level files
        population: true if only to be cleaned at the population level, must include:
            individual files (list of individual-level file paths to be concatenated)
    """    
    # make output directory
    os.system(f'mkdir {output_dir}')
    # get final save path
    save_path = f'{output_dir}/{output_file_name}'
    if individual and input_path:
        # read in the bed file into a pandas dataframe
        df = pd.read_csv(input_path, sep='\t', header=None)
        header = ['Chromosome', 'Start', 'End', 'Name', 'Score', 'Strand', 'ThickStart', 'ThickEnd', 'ItemRGB']
        df.columns = header[:len(df.columns)]
        # keep only denisovan/neanderthal introgression data
        if archaic=='Neanderthal':
            filtered = df[df['Name'].str.contains("neaTo")==True]
        elif archaic=='Denisovan':
            filtered = df[df['Name'].str.contains("denTo")==True]
        # add population and ancestry columns
        filtered['Population']=population
        filtered['Ancestry']=ancestry
        # add archaic introgression source
        filtered["IntrogressionSource"] = filtered.Name.apply(lambda x: "Neanderthal" if "neaTo" in x else ("Denisovan" if "denTo" in x else "Other"))
        # add score name
        filtered["ScoreName"] = 'homozygous or heterozygous'
        # add method name
        filtered["Method"] = 'ArgWeaver-D, Hubisz 2020'
        # add data source
        filtered["GenomesSource"] = genomes_source
        # add individual ID
        filtered['ID'] = ID
        # drop columns that are not of interest
        filtered.drop(['Strand', 'ThickStart', 'ThickEnd', 'ItemRGB'], axis=1, inplace=True)
        # save final file
        filtered.to_csv(f'{output_dir}/{output_file_name}', sep='\t', index=False)
        # completed
        print(f'Output .bed file is located at {save_path}', sep='')

    elif population and individual_files:
        # make directory to store bedgraphs
        os.system(f'mkdir {output_dir}/bedgraphs')
        dfs = []
        for file_path in individual_files:
            # read csv file into a DataFrame
            df = pd.read_csv(file_path, sep='\t')
            # append DataFrame to list
            dfs.append(df)
        # concatenate DataFrames in the list along rows
        df = pd.concat(dfs, ignore_index=True)
        # make directory to store bedgraphs
        os.system(f'mkdir {output_dir}/bedgraphs')
        # convert data to bedgraph format
        save_bedgraph(df, output_dir+'/bedgraphs')
        # get ancestries as a list and string separated by a space
        ancestries_list = list(df['Ancestry'].unique())
        # get list of bedgraph files to be merged in the next step, format is full path separated by " \"
        merged_bg_files = " \\".join([f'{output_dir}'+'/bedgraphs/'+file for file in os.listdir(f'{output_dir}/bedgraphs') if file.startswith("merged_")])
        # get list of ancestries to be used in unionbedg (must be in the same order as the input files)
        ancestries_str = re.sub(f'{output_dir}'+'/bedgraphs/'+r'(merged_|\.bg)', '', merged_bg_files)
        ancestries_str = re.sub(r'(|\.bg)', '', ancestries_str)
        ancestries_str = re.sub(re.escape('\\'), "", ancestries_str)
        # create union file at the population level
        os.system(f'{BEDTOOLSDIR}/bedtools unionbedg -i {merged_bg_files} -header -names {ancestries_str} > {output_dir}/bedgraphs/union_merged.bg')
        # delete intermediate files, keeping only merged files at the ancestral level
        os.system(f'cd {output_dir}/bedgraphs; rm -f sorted*.bg')
        # read in generated union file
        argweaverd_final = pd.read_csv(f'{output_dir}/bedgraphs/union_merged.bg', sep='\t')

        # keep only the highest reported score, collapse into one column
        argweaverd_final['highest_score'] = argweaverd_final[ancestries_list].max(axis=1)
        # drop ancestries
        argweaverd_final = argweaverd_final.drop(ancestries_list, axis=1)
        # save as a bedgraph
        argweaverd_final.to_csv(f'{output_dir}/bedgraphs/union_merged_highest_score.bg', sep='\t', index=False)

        # add metadata
        # add score_name and fill with value "LOD"
        argweaverd_final['ScoreName'] = 'homozygous or heterozygous'
        # add IntrogressionSource
        argweaverd_final['IntrogressionSource'] = archaic
        # add method name
        argweaverd_final['Method'] = 'ArgWeaver-D, Hubisz 2020'
        # add genomes source (user-inputted parameter)
        argweaverd_final['GenomesSource'] = genomes_source

        # save final file
        argweaverd_final.to_csv(f'{output_dir}/{output_file_name}', sep='\t', index=False)
        
        # completed
        print(f'Output .bed file is located at {save_path}', sep='')

def yuan_2021_raw_to_bed(output_dir, output_file_name, archaic, input_dir_list=None, genomes_source=None, file_path=None, population=False, individual=False):
    """
    This function takes in a .seg file (such as IntrogressedSeg/EastAsia/CDX_mergeLK.seg) and cleans it into a standardized output
    Args:
        output_dir: output directory to place the cleaned file 
        output_file_name: output file name to place in the output_dir
        archaic: whether or not to save Neanderthal or Denisovan ancestry, can be 'Neanderthal' or 'Denisovan'
        input_dir_list: list of input directories (to be used with individual=True) 
        genomes_source: the project that the genomes are from (usually 1000 Genomes Project (1KG) or Simons Genome Diversity Project (SGDP))
        file_path = input file containing individual=level introgressed regions (to be used with population=True)
        population: true if to be cleaned at the population level, input_path file should already be cleaned at the individual level
        individual: true if only to be cleaned at the individual level
    """
    if input_dir_list!=None:
        dfs=[] #list to store individual population bed files to be concatenated
        # make output directory
        os.system(f'mkdir {output_dir}')
        for input_dir in input_dir_list:
            if individual:
                # loop through each input file, each representing a different population
                for input_file in os.listdir(input_dir):
                    # must end with .seg
                    if input_file.endswith('.seg'):
                        df = pd.read_csv(f'{input_dir}/{input_file}', sep='\t')
                        # get final save path
                        save_path = f'{output_dir}/{output_file_name}'
                        # filter by if they are denisovan/neanderthal population matches
                        # keep only denisovan/neanderthal introgression data
                        if archaic=='Neanderthal':
                            filtered = df[df['BestMatchedPop']=="Altai"]
                        elif archaic=='Denisovan':
                            filtered = df[df['BestMatchedPop']=="Denisova"]
                        else:
                            print('Must be Neanderthal or Denisovan ancestry')
                        # add population, ancestry, method name, genome source, score name
                        filtered['Population']=input_file.split("_")[0]
                        filtered['Ancestry']=filtered['Population'].apply(map_status)
                        filtered["Method"] = 'ArchaicSeeker 2.0, Yuan 2021'
                        if genomes_source:
                            filtered['GenomesSource'] = genomes_source
                        else: # assume original input from archaicseeker2 paper, which uses 1KG Phase 3 genomes
                            filtered['GenomesSource'] = '1KG and SGDP'
                        filtered["ScoreName"] = 'BestMatchedTime'
                        # change 'Contig' column name to 'Chromosome'
                        filtered.rename(columns={'Contig':'Chromosome', 'Start(bp)':'Start', 'End(bp)':'End', 'BestMatchedTime':'Score'}, inplace=True)
                        # append individual population file to overall file
                        dfs.append(filtered)
        final_bed = pd.concat(dfs)
        # reset index of concatenated df
        final_bed.reset_index(drop=True, inplace=True)
        # drop columns that are not of interest
        final_bed.drop(['Start(cM)', 'End(cM)', 'BestMatchedPop'], axis=1, inplace=True)
        # save to bed file
        save_to_bed(final_bed, f'{save_path}')
    
    if population == True and file_path:
        df = pd.read_csv(file_path, sep='\t')
        # make directory to store bedgraphs
        os.system(f'mkdir {output_dir}/bedgraphs')
        # get final save path
        save_path = f'{output_dir}/{output_file_name}'
        # convert data to bedgraph format
        save_bedgraph(df, output_dir+'/bedgraphs')
        # get ancestries as a list and string separated by a space
        ancestries_list = list(df['Ancestry'].unique())
        # get list of bedgraph files to be merged in the next step, format is full path separated by " \"
        merged_bg_files = " \\".join([f'{output_dir}'+'/bedgraphs/'+file for file in os.listdir(f'{output_dir}/bedgraphs') if file.startswith("merged_")])
        # get list of ancestries to be used in unionbedg (must be in the same order as the input files)
        ancestries_str = re.sub(f'{output_dir}'+'/bedgraphs/'+r'(merged_|\.bg)', '', merged_bg_files)
        ancestries_str = re.sub(r'(|\.bg)', '', ancestries_str)
        ancestries_str = re.sub(re.escape('\\'), "", ancestries_str)
        # create union file at the population level
        os.system(f'{BEDTOOLSDIR}/bedtools unionbedg -i {merged_bg_files} -header -names {ancestries_str} > {output_dir}/bedgraphs/union_merged.bg')
        # delete intermediate files, keeping only merged files at the ancestral level
        os.system(f'cd {output_dir}/bedgraphs; rm -f sorted*.bg')
        # read in generated union file
        archaicseeker2_final = pd.read_csv(f'{output_dir}/bedgraphs/union_merged.bg', sep='\t')

        # keep only the highest reported score, collapse into one column
        archaicseeker2_final['highest_score'] = archaicseeker2_final[ancestries_list].max(axis=1)
        # drop ancestries
        archaicseeker2_final = archaicseeker2_final.drop(ancestries_list, axis=1)
        # save as a bedgraph
        archaicseeker2_final.to_csv(f'{output_dir}/bedgraphs/union_merged_highest_score.bg', sep='\t', index=False)

        # add metadata
        # add score_name and fill with value "LOD"
        archaicseeker2_final['ScoreName'] = 'BestMatchedTime'
        # add IntrogressionSource (Neanderthal or Denisovan, ArchaicSeeker can detect these archaic species)
        archaicseeker2_final['IntrogressionSource'] = archaic
        # add method name
        archaicseeker2_final['Method'] = 'ArchaicSeeker 2.0, Yuan 2021'
        # add genomes source (user-inputted parameter)
        if genomes_source:
            archaicseeker2_final['GenomesSource'] = genomes_source
        else: # assume original input from paper
            archaicseeker2_final['GenomesSource'] = '1KG and SGDP'

        # save final file
        archaicseeker2_final.to_csv(f'{output_dir}/{output_file_name}', sep='\t', index=False)
        
        # completed
        print(f'Output .bed file is located at {save_path}', sep='')

def durvasula_2019_raw_to_bed(output_dir, output_file_name, ancestry, input_dir=None, genomes_source=None, file_path=None, population=False, individual=False, threshold=0.892, keep_old_score=False, keep_average_score=False):
    """
    This function takes in a raw archie output (POP.chr*-prediction.bed.gz-trimmed.bed.gz) and cleans it into a standardized output
    Args:
        output_dir: output directory to place the cleaned file 
        output_file_name: output file name to place in the output_dir
        ancestry: specify ancestry (CEU)
        input_dir: input directory containing individual chromosome introgression files
        genomes_source: the project that the genomes are from (usually 1000 Genomes Project (1KG) or Simons Genome Diversity Project (SGDP))
        file_path = input file containing individual=level introgressed regions (to be used with population=True)
        population: true if to be cleaned at the population level, input_path file should already be cleaned at the individual level
        individual: true if only to be cleaned at the individual level
        threshold: min threshold (Archaic probability) to be considered introgressed, default taken from the paper
        keep_old_score: TRUE keeps the original scores as "Score" and "Score_new" is the heterozygous/homozygous score (I changed Archie's scores to be 0.5 or 1, heterozygous or homozygous). This also keeps only regions above the specified threshold.
        keep_average_score: instead of saving the highest score reported, I'll save the average score for a region
    """
    # Custom function to replace values based on homo (1) or heterozygous (0.5)
    def replace_values(entry):
        values = [float(val) for val in entry.split(';')]
        
        if any(val > threshold for val in values):
            if all(val > threshold for val in values):
                return 1
            else:
                return 0.5
        else:
            return 0
        
    # make output directory
    os.system(f'mkdir {output_dir}')
    if individual and input_dir!=None:
        dfs=[] #list to store individual population bed files to be concatenated
        for input_file in os.listdir(input_dir):
            if input_file != '.DS_Store':
                # read in the bed file
                archie = gzip.open(f'{input_dir}/{input_file}')
                print(f'Reading {input_file}')
                # Will read in as a Pandas dataframe for easy data cleaning but the above line is if you want to read into PyRanges
                archie = pd.read_csv(archie, sep='\t', header=None)
                header = ['Chromosome', 'Start', 'End', 'Score', 'NotArchaicProbability', 'ID', 'Haplotype']
                archie.columns = header[:len(archie.columns)] 
                # add left and right haplotypes based on row number
                archie['Haplotype']=np.where(archie.index%2==0, 'Left', 'Right')
                # Group by each window for each person
                archie = archie.groupby(['ID', 'Chromosome', 'Start', 'End', 'Haplotype'])['Score'].apply(lambda x: ';'.join(x.astype(str))).unstack().reset_index()
                # make sure columns are floats
                archie['Left'] = pd.to_numeric(archie['Left'], errors='coerce')
                archie['Right'] = pd.to_numeric(archie['Right'], errors='coerce')
                # keep only scores greater than a threshold (default 0.5)
                archie_filt = archie[(archie['Left'] >= threshold) | (archie['Right'] >= threshold)]
                # convert back to strings
                archie_filt['Left'] = archie_filt['Left'].apply(str)
                archie_filt['Right'] = archie_filt['Right'].apply(str)
                # Combine scores into semicolon-separated format
                archie_filt['Score'] = archie_filt['Left'].fillna('') + ';' + archie_filt['Right'].fillna('')
                # Split the 'Score' column into a list
                archie_filt['Score'] = archie_filt['Score'].str.split(';')
                # Use explode to create separate rows
                archie_filt_expanded = archie_filt.explode('Score', ignore_index=True)
                # Convert the Values column back to float
                archie_filt_expanded['Score'] = archie_filt_expanded['Score'].astype(float)
                # Drop uncombined columns
                archie_filt = archie_filt_expanded.drop(['Left', 'Right'], axis=1)
                # add population
                archie_filt['Ancestry']=ancestry
                # add method name
                archie_filt['Method']='ArchIE'
                # add GenomesSource
                archie_filt['GenomesSource']=genomes_source
                # add ScoreName
                archie_filt['ScoreName']='ArchaicProbability'
                # remove "name" attribute
                archie_filt.columns.name = None
                # Apply the custom function to the specified column if keep_old_score = False (default)
                if keep_old_score == False:
                    archie_filt['Score_new'] = archie_filt['Score'].apply(replace_values)
                    # Rename score column
                    archie_filt = archie_filt.drop(columns=['Score'])
                    archie_filt.rename(columns={'Score_new': 'Score'}, inplace=True)
                    # Add ScoreName
                    archie_filt['ScoreName'] = f'Thresholded archaic probability, >{threshold} considered introgressed in the raw files'
                    # subset df so we only have scores above the threshold
                    archie_filt = archie_filt[archie_filt['Score'] > threshold]
                else:
                    # subset df so we only have scores above the threshold
                    archie_filt = archie_filt[archie_filt['Score'] > threshold]
                    # Add ScoreName
                    archie_filt['ScoreName_New'] = f'0.5 for heterzygous, 1 for homozygous. Thresholded archaic probability, >{threshold} considered introgressed in the raw files'
                    archie_filt['ScoreName'] = 'Introgression Probability'
                # NOTE
                # Archie uses a different numbering system than UCSC genome browser, need to subtract 1 from the end for all files
                # yes that means we only have 49999 base window sizes but then we avoid double counting a base twice
                archie_filt['End'] = archie_filt['End'] - 1
                # append individual population file to overall file
                dfs.append(archie_filt)
        # concatenate the dataframes
        final_bed = pd.concat(dfs)
        # reset index of concatenated df
        final_bed.reset_index(drop=True, inplace=True)
        # save to bed file
        # get final save path
        save_path = f'{output_dir}/{output_file_name}'
        final_bed.to_csv(f'{save_path}', sep='\t', index=False, header=True)
        #save_to_bed(final_bed, f'{save_path}')
            
    if population == True and file_path:
        df = pd.read_csv(file_path, sep='\t')
        # make directory to store bedgraphs
        os.system(f'mkdir {output_dir}/bedgraphs')
        # get final save path
        save_path = f'{output_dir}/{output_file_name}'
        # convert data to bedgraph format
        save_bedgraph(df, output_dir+'/bedgraphs')
        # get ancestries as a list and string separated by a space
        ancestries_list = list(df['Ancestry'].unique())
        # if more than one ancestry, merge individual .bg file
        # get list of bedgraph files to be merged in the next step, format is full path separated by " \"
        merged_bg_files = " \\".join([f'{output_dir}'+'/bedgraphs/'+file for file in os.listdir(f'{output_dir}/bedgraphs') if file.startswith("merged_")])
        # get list of ancestries to be used in unionbedg (must be in the same order as the input files)
        ancestries_str = re.sub(f'{output_dir}'+'/bedgraphs/'+r'(merged_|\.bg)', '', merged_bg_files)
        ancestries_str = re.sub(r'(|\.bg)', '', ancestries_str)
        ancestries_str = re.sub(re.escape('\\'), "", ancestries_str)
        if len(ancestries_list)>1: # if more than one ancestry group, merge them
            # create union file at the population level
            os.system(f'{BEDTOOLSDIR}/bedtools unionbedg -i {merged_bg_files} -header -names {ancestries_str} > {output_dir}/bedgraphs/union_merged.bg')
        else: # rename merged population file to union_merged.bg
            os.system(f'cd {output_dir}/bedgraphs; cp merged_{ancestries_str}.bg union_merged.bg')            
        # delete intermediate files, keeping only merged files at the ancestral level
        os.system(f'cd {output_dir}/bedgraphs; rm -f sorted*.bg')
        [os.remove(os.path.join(f'{output_dir}/bedgraphs', file)) for file in os.listdir(f'{output_dir}/bedgraphs') if any(file.startswith(ancestry) for ancestry in ancestries_list)]
        # read in generated union file
        archie_final = pd.read_csv(f'{output_dir}/bedgraphs/union_merged.bg', sep='\t', names=['chrom', 'start', 'end', f'{ancestry}'])
        # save this version with the proper header
        archie_final.to_csv(f'{output_dir}/bedgraphs/union_merged.bg', sep='\t', index=False)
        # keep only the highest reported score, collapse into one column
        archie_final['highest_score'] = archie_final[ancestries_list].max(axis=1)
        # drop ancestries
        archie_final = archie_final.drop(ancestries_list, axis=1)
        # save as a bedgraph
        archie_final.to_csv(f'{output_dir}/bedgraphs/union_merged_highest_score.bg', sep='\t', index=False, header=False)

        # add metadata
        archie_final['ScoreName'] = f'{threshold} for heterzygous, 1 for homozygous. Thresholded archaic probability, >{threshold} considered introgressed in the raw files'
        # add IntrogressionSource, Neanderthal for ArchIE
        archie_final['IntrogressionSource'] = 'Neanderthal'
        # add method name
        archie_final['Method'] = 'ArchIE, Durvasula 2019'
        # add genomes source (user-inputted parameter)
        if genomes_source:
            archie_final['GenomesSource'] = genomes_source
        else: # assume original input from archie paper, which uses 1KG CEU genomes
            archie_final['GenomesSource'] = '1KG CEU'
        # save final file
        archie_final.to_csv(f'{output_dir}/{output_file_name}', sep='\t', index=False)

        # completed
        print(f'Output .bed file is located at {save_path}', sep='')

def steinruecken_2018_raw_to_bed(output_dir, output_file_name, input_tar_gz=None, input_dir=None, prob_dir=None, genomes_source=None, file_path=None, population=False, individual=False):
    """
    This function takes in an input directory containing .tar.gz files for each population/chromosome and cleans it to a standardized output.
    HOWEVER, I highly recommend running this through parallel bash scripts. See submit_cleansteinrucken_parallel.sh
    Thus, this function is meant to be run for one .tar.gz file.
    Args:
        output_dir: output directory to place the cleaned file, will contain cleaned population/chromosome-level files if individual=True
        output_file_name: output file name to place in the output_dir
        input_dir_list: list of input directories (to be used with individual=True)
        prob_dir: posterior probabilities directory
        genomes_source: the project that the genomes are from (usually 1000 Genomes Project (1KG) or Simons Genome Diversity Project (SGDP))
        file_path = input file containing individual=level introgressed regions (to be used with population=True)
        population: true if to be cleaned at the population level, input_path file should already be cleaned at the individual level
        individual: true if only to be cleaned at the individual level, then must include:
            input_dir, prob_dir, output_file_name (final path of file), input_tar_gz
        population: true if only to be cleaned at the population level, then must include:
            file_path (output of individual=True), output_dir, output_file_name (final path of file)
    """
    def collapse_columns(df,col):
        cols = df.columns.values.tolist()
        cols.remove(col)
        df = df.groupby(cols)[col].apply(';'.join).reset_index()
        return df
    if individual and input_dir:
        # get ancestries list from input file names, keeping unique ancestries
        ancestries_list = list(set([file.split('_')[0] for file in os.listdir(input_dir) if file.endswith('_beds.tar.gz')]))
        # overall df
        steinruecken = pd.DataFrame() 
        # make output directory
        os.system(f'mkdir {output_dir}')
        #loop over each file
        for i in ancestries_list:
            # get the tar.gz file passed in
            j = input_dir + '/' + input_tar_gz
        # loop over ancestries
        if j.startswith(i):
            f = input_dir + j
            PROB_DIR = prob_dir + '/' + j.replace('_beds', '')
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
                position_probs.columns=['Start','Score']
                df_temp.columns=['Chromosome','Start', 'End', 'Ancestry', 'ID']
                df_temp[['Chromosome', 'Ancestry', 'ID']] = df_temp[['Chromosome', 'Ancestry', 'ID']].astype(str)
                # merge the data only if the positions appear in both
                final = df_temp.merge(position_probs, on='Start')
                steinruecken = pd.concat([steinruecken, final])

                OUTPUT_FILE = f'{output_dir}/{j}_frag.bed'
                steinruecken.rename(columns={0:'Chr',1:'Start',2:'End', 3:'Ancestry', 4:'ID', 5:'Score'}, inplace=True)
                steinruecken.to_csv(OUTPUT_FILE, sep='\t', index=False)
                # Merge all the files and save as a single bed file
                tracts = os.listdir(output_dir)
                tracts = [i for i in tracts if i.endswith('frag.bed')]
                groups = list(set([re.sub(r'([A-Z]+)_.*',r'\1',i) for i in tracts]))
                # initialize empty dataframe
                steinruecken = pd.DataFrame()
                # concatenate the files
                for i in groups:
                    for j in tracts:
                        if j.startswith(i):
                            df_temp = pd.read_csv(output_dir + '/' + j, sep='\t', header=0)
                            steinruecken = pd.concat([steinruecken, df_temp])

                # drop duplicates
                steinruecken = steinruecken.drop_duplicates()
                # save to csv
                steinruecken.to_csv(output_file_name, sep='\t', index=False)
    
    if population == True and file_path:
        df = pd.read_csv(file_path, sep='\t')
        # make directory to store bedgraphs
        os.system(f'mkdir {output_dir}/')
        os.system(f'mkdir {output_dir}/bedgraphs')
        # get final save path
        save_path = f'{output_dir}/{output_file_name}'
        # convert data to bedgraph format
        save_bedgraph(df, output_dir+'/bedgraphs')
        # get ancestries as a list and string separated by a space
        ancestries_list = list(df['Ancestry'].unique())
        # get list of bedgraph files to be merged in the next step, format is full path separated by " \"
        merged_bg_files = " \\".join([f'{output_dir}'+'/bedgraphs/'+file for file in os.listdir(f'{output_dir}/bedgraphs') if file.startswith("merged_")])
        # get list of ancestries to be used in unionbedg (must be in the same order as the input files)
        ancestries_str = re.sub(f'{output_dir}'+'/bedgraphs/'+r'(merged_|\.bg)', '', merged_bg_files)
        ancestries_str = re.sub(r'(|\.bg)', '', ancestries_str)
        ancestries_str = re.sub(re.escape('\\'), "", ancestries_str)
        # create union file at the population level
        os.system(f'{BEDTOOLSDIR}/bedtools unionbedg -i {merged_bg_files} -header -names {ancestries_str} > {output_dir}/bedgraphs/union_merged.bg')
        # delete intermediate files, keeping only merged files at the ancestral level
        os.system(f'cd {output_dir}/bedgraphs; rm -f sorted*.bg')
        # read in generated union file
        dicaladmix_final = pd.read_csv(f'{output_dir}/bedgraphs/union_merged.bg', sep='\t')
        # remove 'chr' prefix
        dicaladmix_final['chrom'] = dicaladmix_final['chrom'].str.replace('chr', '')

        # keep only the highest reported score, collapse into one column
        dicaladmix_final['highest_score'] = dicaladmix_final[ancestries_list].max(axis=1)
        # drop ancestries
        dicaladmix_final = dicaladmix_final.drop(ancestries_list, axis=1)
        # save as a bedgraph
        dicaladmix_final.to_csv(f'{output_dir}/bedgraphs/union_merged_highest_score.bg', sep='\t', index=False)

        # add metadata
        # add score_name and fill with value "LOD"
        dicaladmix_final['ScoreName'] = 'Posterior probabilities'
        # add IntrogressionSource (Neanderthal)
        dicaladmix_final['IntrogressionSource'] = 'Neanderthal'
        # add method name
        dicaladmix_final['Method'] = 'DICAL-ADMIX, Steinruecken 2018'
        # add genomes source (user-inputted parameter)
        if genomes_source:
            dicaladmix_final['GenomesSource'] = genomes_source
        else: # assume original input from paper
            dicaladmix_final['GenomesSource'] = '1000 Genomes Phase 1'
        # remove 'chr' prefix
        dicaladmix_final['chrom'] = dicaladmix_final['chrom'].str.replace('chr', '')

        # save final file
        dicaladmix_final.to_csv(f'{output_dir}/{output_file_name}', sep='\t', index=False)
        # completed
        print(f'Output .bed file is located at {save_path}', sep='')

def skov_2020_raw_to_bed(input_file, output_dir, output_file_name, archaic, population=False, ancestry='Icelandic', liftover_path=None, genomes_source='Icelandic Genomes', snp_level=False):
    """
    This function takes in a raw skov 2020 output (41586_2020_2225_MOESM3_ESM.txt from Supplementary material) and cleans it into a standardized output
    Args:
        input_file: input file containing introgression data (41586_2020_2225_MOESM3_ESM.txt)
        liftover_path: include liftover path to run liftover to convert hg38 to hg19. the liftover directory should contain the hg39 to hg19 chain.
        output_dir: output directory to place the cleaned file 
        archaic: extracting Neanderthal or Denisovan introgression
        output_file_name: output file name to place in the output_dir
        genomes_source: the project that the genomes are from (usually 1000 Genomes Project (1KG) or Simons Genome Diversity Project (SGDP))
        file_path = input file containing individual=level introgressed regions (to be used with population=True)
        population: true if to be cleaned at the population level, input_path file should already be cleaned at the individual level
        ancestry: ancestry of the individual represented, default is 'Icelandic'
        snp_level: clean the snps file. If true, the input file should be 41586_2020_2225_MOESM4_ESM.txt
    """
    if snp_level:
        skov_2020_snps = pd.read_csv(f'{input_file}', 
                        sep='\t', 
                        low_memory=False)
        skov_2020_snps['start'] = skov_2020_snps['pos'] - 1
        # keep only DAV snps
        skov_2020_snps = skov_2020_snps[(skov_2020_snps['snptype'] == 'DAVsnp') | (skov_2020_snps['snptype'] == 'linkedDAVsnp')]
        # rename pos to end
        skov_2020_snps.rename(columns={'pos': 'end'}, inplace=True)
                # filter by archiac segments of interest
        if archaic == 'Neanderthal':
            skov_2020_snps_filt = skov_2020_snps[(skov_2020_snps['foundin'].str.contains('Vindija')) | (skov_2020_snps['foundin'].str.contains('Altai'))]
        elif archaic == 'Denisovan':
            skov_2020_snps_filt = skov_2020_snps[(skov_2020_snps['foundin'] == 'Denisova')]
        # reorder extra columns
        skov_2020_snps_filt=skov_2020_snps_filt[['chrom', 'start', "end", "ref", "alt", "derived", "freqinPOP", "freqinArchaicsegments", "foundin", "snptype", "DAVmarker", "DAVLD", "nonDAVmarker", "nonDAVLD", "maxLD"]]
        skov_2020_snps_filt.to_csv(f'{output_dir}/snps_for_liftover.bed', sep='\t', header=False, index=False)        
        # run liftover
        os.system(f'cd {liftover_path}; ./liftOver {output_dir}/snps_for_liftover.bed {liftover_path}/hg38ToHg19.over.chain.gz {output_dir}/{output_file_name} {output_dir}/unlifted_snps_for_liftover.bed -bedPlus=3')
        # add a header
        os.system(f"sed -i $'1ichrom\tstart\tend\tref\talt\tderived\tfreqinPOP\tfreqinArchaicsegments\tfoundin\tsnptype\tDAVmarker\tDAVLD\tnonDAVmarker\tnonDAVLD\tmaxLD' {output_dir}/{output_file_name}")
        
    elif population:
        # read in full file
        skov_2020 = pd.read_csv(f'{input_file}', sep='\t', low_memory=False)
        skov_2020.columns = ['chrom', 'start', 'end', 'type', 'length', 'DAVsnps', 'altai', 'denisova', 'vindija', 'totalsnps', 'archaic', 'called', 'freq', 'privateAltai', 'privateDenisova', 'privateVindija', 'mosaicstatus', 'snp_list']

        if liftover_path:
            # write to output
            skov_2020.to_csv(f'{output_dir}/skov_2020_for_liftover.txt', sep='\t', index=False, header=False)
            # liftover to hg19
            os.system(f'./{liftover_path}/liftover {output_dir}/skov_2020_for_liftover.txt hg38ToHg19.over.chain.gz {output_dir}/skov_2020_lifted.bed {output_dir}/unlifted.bed  -bedPlus=3')

            # open the hg19-lifted variant file
            skov_2020 = pd.read_csv(f'{output_dir}/skov_2020_lifted.bed', sep='\t', low_memory=False)
            skov_2020.columns = ['chrom', 'start', 'end', 'type', 'length', 'DAVsnps', 'altai', 'denisova', 'vindija', 'totalsnps', 'archaic', 'called', 'freq', 'privateAltai', 'privateDenisova', 'privateVindija', 'mosaicstatus', 'snp_list']
            # remove 'chr' prefix
            skov_2020['chrom'] = skov_2020['chrom'].str.replace('chr', '')

        # filter by archiac segments of interest
        if archaic == 'Neanderthal':
            skov_2020_filtered = skov_2020[(skov_2020['archaic'] == 'Vindija') | (skov_2020['archaic'] == 'Altai')]
        if archaic == 'Denisovan':
            skov_2020_filtered = skov_2020[(skov_2020['archaic'] == 'Denisova')]

        # add IntrogressionSource 
        skov_2020_filtered['IntrogressionSource'] = archaic
        # add GenomesSource 
        skov_2020_filtered['Ancestry'] = 'Icelandic'
        # drop irrelevant columns
        skov_2020_filtered.drop(['type', 'length', 'privateAltai', 'privateDenisova', 'privateVindija', 'mosaicstatus', 'DAVsnps', 'altai', 'denisova', 'vindija', 'totalsnps'], axis=1, inplace=True)
        # rename columns
        skov_2020_filtered.rename(columns={'start': 'Start'}, inplace=True)
        skov_2020_filtered.rename(columns={'end': 'End'}, inplace=True)
        skov_2020_filtered.rename(columns={'chrom': 'Chromosome'}, inplace=True)
        skov_2020_filtered['Score'] = 1
        # remove 'chr' prefix
        if skov_2020_filtered['Chromosome'].str.contains('chr').any():
            skov_2020_filtered['Chromosome'].replace('[chr]', '', inplace=True, regex=True)

        # make directory to store bedgraphs
        os.system(f'mkdir {output_dir}/')
        os.system(f'mkdir {output_dir}/bedgraphs')
        # get final save path
        save_path = f'{output_dir}/{output_file_name}'
        # convert data to bedgraph format
        save_bedgraph(skov_2020_filtered, output_dir+'/bedgraphs')
        # get ancestries as a list and string separated by a space
        ancestries_list = list(skov_2020_filtered['Ancestry'].unique()) # should be just Icelandic
        # if more than one ancestry, merge individual .bg file
        # get list of bedgraph files to be merged in the next step, format is full path separated by " \"
        merged_bg_files = " \\".join([f'{output_dir}'+'/bedgraphs/'+file for file in os.listdir(f'{output_dir}/bedgraphs') if file.startswith("merged_")])
        # get list of ancestries to be used in unionbedg (must be in the same order as the input files)
        ancestries_str = re.sub(f'{output_dir}'+'/bedgraphs/'+r'(merged_|\.bg)', '', merged_bg_files)
        ancestries_str = re.sub(r'(|\.bg)', '', ancestries_str)
        ancestries_str = re.sub(re.escape('\\'), "", ancestries_str)
        if len(ancestries_list)>1: # if more than one ancestry group, merge them
            # create union file at the population level
            os.system(f'{BEDTOOLSDIR}/bedtools unionbedg -i {merged_bg_files} -header -names {ancestries_str} > {output_dir}/bedgraphs/union_merged.bg')
        else: # rename merged population file to union_merged.bg
            os.system(f'cd {output_dir}/bedgraphs; cp merged_{ancestries_str}.bg union_merged.bg')

        # delete intermediate files, keeping only merged files at the ancestral level
        os.system(f'cd {output_dir}/bedgraphs; rm -f sorted*.bg')
        [os.remove(os.path.join(f'{output_dir}/bedgraphs', file)) for file in os.listdir(f'{output_dir}/bedgraphs') if any(file.startswith(ancestry) for ancestry in ancestries_list)]
        # read in generated union file
        skov2020_final = pd.read_csv(f'{output_dir}/bedgraphs/union_merged.bg', sep='\t', names=['chrom', 'start', 'end', f'{ancestry}'])

        # keep only the highest reported score, collapse into one column
        skov2020_final['highest_score'] = skov2020_final[ancestries_list].max(axis=1)
        # drop ancestries
        skov2020_final = skov2020_final.drop(ancestries_list, axis=1)
        # save as a bedgraph
        skov2020_final.to_csv(f'{output_dir}/bedgraphs/union_merged_highest_score.bg', sep='\t', index=False)

        # add metadata
        skov2020_final['ScoreName'] = 'Default of 1 since all fragments have a prob > 0.9'
        # add IntrogressionSource, Neanderthal or Denisovan
        skov2020_final['IntrogressionSource'] = archaic
        # add method name
        skov2020_final['Method'] = 'HMM, Skov 2020'
        # add genomes source (user-inputted parameter)
        if genomes_source:
            skov2020_final['GenomesSource'] = genomes_source
        else: # assume original input from archie paper, which uses 1KG CEU genomes
            skov2020_final['GenomesSource'] = 'Icelandic Genomes'

        # change ancestry back
        skov2020_final.rename(columns={'highest_score': 'Icelandic'}, inplace=True)
        # save final file
        skov2020_final.to_csv(save_path, sep='\t', index=False)
        # completed
        print(f'Output .bed file is located at {save_path}', sep='')

def schaefer_2021_raw_to_bed(input_file, output_dir, output_file_name, archaic, population=None, individual=None, genomes_source='SGDP'):
    """
    This method reads the raw output from schaefer_2021 (sarge) and cleans it to .bed format with overlapping regions collapsed into one entry, either at the individual or population-levels.
    Args:
        input_file: path to raw SARGE file,
        output_dir: output directory to place the cleaned file 
        output_file_name: output file name to place in the output_dir
        archaic: which archaic species to include ('Neanderthal' or 'Denisovan')
        genomes_source: the source that the genomes are from (usually 1000 Genomes Project (1KG) or Simons Genome Diversity Project (SGDP))
        population: true if to be cleaned at the population level, input_path file should already be cleaned at the individual level
        individual: true if only to be cleaned at the individual level
    """
    df = pd.read_csv(input_file, sep='\t')
    # make directory to store bedgraphs
    os.system(f'mkdir {output_dir}')
    os.system(f'mkdir {output_dir}/bedgraphs')
    # get final save path
    save_path = f'{output_dir}/{output_file_name}'
    # keep only fragments of interest
    if archaic=='Neanderthal':
        df = df[(df['type']=='NEA')]
    elif archaic=='Denisovan':
        df = df[(df['type']=='DEN')]
    # drop unneeded columns
    df.drop(['hap', 'mean_tmrca', 'len'], axis=1, inplace=True)
    # add columns of interest
    df['Method'] = 'SARGE'
    df['GenomesSource'] = 'SGDP'
    df['IntrogressionSource'] = 'Neanderthal'
    df.rename(columns={'pop': 'Ancestry'}, inplace=True)
    df.rename(columns={'#chrom': 'Chromosome'}, inplace=True)
    df.rename(columns={'start': 'Start'}, inplace=True)
    df.rename(columns={'end': 'End'}, inplace=True)
    df.rename(columns={'score': 'Score'}, inplace=True)
    if population == True:
        # convert data to bedgraph format
        save_bedgraph(df, output_dir+'/bedgraphs')
        # get ancestries as a list and string separated by a space
        ancestries_list = list(df['Ancestry'].unique())
        # get list of bedgraph files to be merged in the next step, format is full path separated by " \"
        merged_bg_files = " \\".join([f'{output_dir}'+'/bedgraphs/'+file for file in os.listdir(f'{output_dir}/bedgraphs') if file.startswith("merged_")])
        # get list of ancestries to be used in unionbedg (must be in the same order as the input files)
        ancestries_str = re.sub(f'{output_dir}'+'/bedgraphs/'+r'(merged_|\.bg)', '', merged_bg_files)
        ancestries_str = re.sub(r'(|\.bg)', '', ancestries_str)
        ancestries_str = re.sub(re.escape('\\'), "", ancestries_str)
        # create union file at the population level
        os.system(f'{BEDTOOLSDIR}/bedtools unionbedg -i {merged_bg_files} -header -names {ancestries_str} > {output_dir}/bedgraphs/union_merged.bg')
        # delete intermediate files, keeping only merged files at the ancestral level
        os.system(f'cd {output_dir}/bedgraphs; rm -f sorted*.bg')
        # read in generated union file
        sarge_final = pd.read_csv(f'{output_dir}/bedgraphs/union_merged.bg', sep='\t')

        # keep only the highest reported score, collapse into one column
        sarge_final['highest_score'] = sarge_final[ancestries_list].max(axis=1)
        # drop ancestries
        sarge_final = sarge_final.drop(ancestries_list, axis=1) #-- keeping this line removes info on which ancestry the introgressed fragment was found in
        # save as a bedgraph
        sarge_final.to_csv(f'{output_dir}/bedgraphs/union_merged_highest_score.bg', sep='\t', index=False)

        # add metadata
        # add score_name
        sarge_final['ScoreName'] = 'Similarity scores'
        # add IntrogressionSource
        sarge_final['IntrogressionSource'] = archaic
        # add method name
        sarge_final['Method'] = 'SARGE, Schaefer 2020'
        # add genomes source (user-inputted parameter)
        sarge_final['GenomesSource'] = genomes_source

        # save final file
        sarge_final.to_csv(f'{output_dir}/{output_file_name}', sep='\t', index=False)
        
        # completed
        print(f'Output .bed file is located at {save_path}', sep='')

def browning_2018_raw_to_bed(input_dir, output_dir, output_file_name, archaic, genomes_source='1KG Phase 3 version 5a', population=True, snp_level=False):
    """
    This method reads the raw output from browning_2018 (sprime) and cleans it to .bed format with overlapping regions collapsed into one entry at the superpopulation-levels.
    Args:
        input_dir: path to raw sprime directory, containing POP_sprime_results.tar.gz
        output_dir: output directory to place the cleaned file 
        output_file_name: output file name to place in the output_dir
        archaic: which archaic species to include ('Neanderthal' or 'Denisovan')
        genomes_source: the source that the genomes are from (usually 1000 Genomes Project (1KG) or Simons Genome Diversity Project (SGDP))
        population: true if to be cleaned at the population level (this is the only option)
        snp_level: true if you want to save the variants, not the fragments
    """
    def create_bed_format(df):
        # add 'chr' prefix to the chromosome names
        # convert to int to remove decimals, and then to string to add 'chr'
        df['Chr'] = df['Chr'].astype(int).astype(str).replace('^', 'chr', regex=True)
        # rename 'POS' column to 'End' and change datatype to int
        #df.rename(columns={'POS':'End'}, inplace=True)
        df['End'] = df['End'].astype(int)
        # create 'Start' column by subtracting 1 from 'End' positions
        df.insert(1, 'Start', df['End']-1)
        
        return df

    def sort_chromosome(df):
        # define sorting order for 'Chr'
        order = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
        # convert 'Chr' to categorical
        df['Chr'] = pd.Categorical(df['Chr'], categories=order, ordered=True)
        # sort by 'Chr' and 'Start'
        df.sort_values(by=['Chr', 'Start'], inplace=True)

        return df
    
    def create_fragment_data(df):
        frags = pd.DataFrame(columns = ['Chromosome', 'Start', 'End', 'Population', 'Score'])
        for chrom in df['Chr'].unique().tolist():
            print(chrom)
            for pop in df['Population'].unique().tolist():
                subset = df[(df['Chr'] == chrom) & (df['Population'] == pop)]
                for segment in subset['Segment_ID'].unique():
                    subset_segment = subset[(subset['Segment_ID'] == segment)]
                    if subset_segment.isna().any().any() == False:
                        start = subset_segment['Start'].min()
                        end = subset_segment['End'].max()
                        frags.loc[len(frags)]=[chrom, start, end, pop, subset_segment['Score'].unique()[0]]
        frags['Method']='Sprime, Browning 2018'
        # add GenomesSource
        frags['GenomesSource']='1KG Phase 3 version 5a'
        # add ScoreName
        frags['ScoreName']='S'
        return frags
    
    # make directory to store data
    os.system(f'mkdir {output_dir}')
    # list files in input_dir
    files = os.listdir(input_dir)
    # parse and merge individual files for each population, including papuans
    dfs = []
    for filename in files:
        if filename.endswith('_sprime_results.tar.gz'):
            # extract population name from the filename
            pop = filename.split('_')[0]
            # get path name
            f = input_dir + '/' + filename
            # open tarfile
            tar = tarfile.open(f, 'r:gz')
            # get chromosomal files
            chr_files = [f.name for f in tar.getmembers()]
            # extract chromosomal data
            for chr_file in chr_files:
                contents = tar.extractfile(chr_file).read()
                f = pd.read_csv(io.BytesIO(contents), sep='\t', header=None, skiprows=1)
                # insert population name
                f.insert(5, 'Population', pop)
                # rename chromosome column
                f = f.rename(columns={f.columns[0]:'Chr', 
                                f.columns[1]:'End', 
                                f.columns[2]:'RSID', 
                                f.columns[3]:'Alt',
                                f.columns[4]:'Ref',
                                f.columns[6]:'Segment_ID', 
                                f.columns[8]:'Score',
                                f.columns[9]:'NMatch',
                                f.columns[10]:'DMatch'})
                # append dataframe to list
                dfs.append(f)
    
    # concatenate dataframes
    snps = pd.concat(dfs, ignore_index=True)
    # remove rows containing empty Chr values
    snps = snps[~snps['Chr'].isna()]
    # sort rows by locus
    snps = snps.sort_values(by=['Chr','End'])
    # wrangle dataframe into bed format
    snps = create_bed_format(snps)
    # separate by archaic hominim
    if archaic=='Neanderthal':
        snps_final = snps[snps['NMatch'] == 'match']
        snps_final = snps_final.drop(['DMatch', 'NMatch'], axis=1)
    if archaic=='Denisovan':
        snps_final = snps[snps['DMatch'] == 'match']
        snps_final = snps_final.drop(['DMatch', 'NMatch'], axis=1)
    # collapse 'Population' column and sort
    cols = ['Chr', 'Start', 'End', 'RSID', 'Ref', 'Alt', 'Segment_ID', 'Score']
    # combine repeated entries
    snps_final_u = sort_chromosome(snps_final)
    # save snps
    if snp_level==True:
        snps_final_u['Chr'] = snps_final_u['Chr'].str.replace('chr', '')
        snps_final_u['Ancestry'] = snps_final_u['Population'].apply(map_status)
        snps_final_u.to_csv(f'{output_dir}/{output_file_name}', sep='\t', index=False)
    else:
        # make dir for bedgraphs
        os.system(f'mkdir {output_dir}/bedgraphs')
        # generate final fragment files
        frags = create_fragment_data(snps_final_u)
        # convert population to ancestry info
        frags['Ancestry'] = frags['Population'].apply(map_status)
        df = frags.drop(columns=['Population'])
        # convert data to bedgraph format
        save_bedgraph(df, output_dir+'/bedgraphs')
        # get ancestries as a list and string separated by a space
        ancestries_list = list(df['Ancestry'].unique())
        # get list of bedgraph files to be merged in the next step, format is full path separated by " \"
        merged_bg_files = " \\".join([f'{output_dir}'+'/bedgraphs/'+file for file in os.listdir(f'{output_dir}/bedgraphs') if file.startswith("merged_")])
        # get list of ancestries to be used in unionbedg (must be in the same order as the input files)
        ancestries_str = re.sub(f'{output_dir}'+'/bedgraphs/'+r'(merged_|\.bg)', '', merged_bg_files)
        ancestries_str = re.sub(r'(|\.bg)', '', ancestries_str)
        ancestries_str = re.sub(re.escape('\\'), "", ancestries_str)
        # create union file at the population level
        os.system(f'{BEDTOOLSDIR}/bedtools unionbedg -i {merged_bg_files} -header -names {ancestries_str} > {output_dir}/bedgraphs/union_merged.bg')
        # delete intermediate files, keeping only merged files at the ancestral level
        os.system(f'cd {output_dir}/bedgraphs; rm -f sorted*.bg')
        # read in generated union file
        sprime_final = pd.read_csv(f'{output_dir}/bedgraphs/union_merged.bg', sep='\t')
        # keep only the highest reported score, collapse into one column
        sprime_final['highest_score'] = sprime_final[ancestries_list].max(axis=1)
        # drop ancestries
        sprime_final = sprime_final.drop(ancestries_list, axis=1)
        # save as a bedgraph
        sprime_final['chrom'] = sprime_final['chrom'].str.replace('chr', '')
        sprime_final.to_csv(f'{output_dir}/bedgraphs/union_merged_highest_score.bg', sep='\t', index=False)
        # add metadata
        # add IntrogressionSource
        sprime_final['IntrogressionSource'] = archaic
        # add method name
        sprime_final['Method']='Sprime, Browning 2018'
        # add GenomesSource
        sprime_final['GenomesSource']=genomes_source
        # add ScoreName
        sprime_final['ScoreName']='S'
    
        # save final file
        sprime_final.to_csv(f'{output_dir}/{output_file_name}', sep='\t', index=False)
    
    # completed
    print(f'Output .bed file is located at {output_dir}/{output_file_name}', sep='')

def vernot_2016_raw_to_bed(raw_introgressed_haplotypes_dir, introgressed_haplotypes_file, output_dir, output_file_name, archaic, genomes_source='1KG and Melanesian samples', individual_level=False):
    """
        This method reads the raw output from vernot 2016 (s star) and cleans it to .bed format with overlapping regions collapsed into one entry at the superpopulation-levels.
        Args:
            raw_introgressed_haplotypes_dir: the path to the dir containing LL.callset$POP.mr_0.99. files (downloaded from untarred tarball (introgressed_haplotypes.tar.gz), https://drive.google.com/drive/folders/0B9Pc7_zItMCVWUp6bWtXc2xJVkk?resourcekey=0-Cj8G4QYndXQLVIGPoWKUjQ)
            introgressed_haplotypes_file is the concatenated files from untarred tarball (introgressed_haplotypes.tar.gz) downloaded from https://drive.google.com/drive/folders/0B9Pc7_zItMCVWUp6bWtXc2xJVkk?resourcekey=0-Cj8G4QYndXQLVIGPoWKUjQ. the introgressed_haplotypes_file is concatenated data from each population's LL.callset"$i".mr_0.99.den_calls_by_hap.bed.merged.by_chr.bed file. Comes from the introgression_methods/bin/vernot_2016_introgressed_haplotypes.sh script as neanderthal_introgressed_haplotypes_individual.bed, containing individual-level haplotype information
            output_dir: output directory to place the cleaned file 
            output_file_name: output file name to place in the output_dir
            archaic: which archaic species to include ('Neanderthal' or 'Denisovan')
            genomes_source: the source that the genomes are from (usually 1000 Genomes Project (1KG) or Simons Genome Diversity Project (SGDP))
            individual_level: True if you're trying to generate an individual-level bed file from the raw directory, to be used with raw_introgressed_haplotypes_dir! The output of this, introgressed_haplotypes_file, should be used to generate the population-level file
    """
    if individual_level == True:
        
        # make list for pop-subset dfs
        dfs=[]
        
        # make temp dir
        os.system(f'mkdir {output_dir}/tmp/')
        
        # for each pop, read in the positions only and positions + raw scores files
        for pop in ['EAS', 'EUR', 'PNG', 'SAS']:
            # positions only (merged)
            haps_no_score = pd.read_csv(f'{raw_introgressed_haplotypes_dir}/LL.callset{pop}.mr_0.99.neand_calls_by_hap.bed.merged.by_chr.bed',
                                sep='\t', low_memory=False, names=['Chromosome', 'Start', 'Stop', 'HaplotypeID', 'ID'])
            # add chr prefix and population info
            haps_no_score['Chromosome'] = 'chr' + haps_no_score['Chromosome'].astype(str)
            
            # read in unmerged positions with scores
            haps_score = pd.read_csv(f'{raw_introgressed_haplotypes_dir}/LL.callset{pop}.mr_0.99.neand_calls_by_hap.bed',
                                sep='\t', low_memory=False, skiprows=1, names=['HaplotypeID', 'Start', 'Stop', 'Chromosome', 'ID', 'Haplotype', 'logit.nean', 'logit.den', 'post.p.null', 'ProbabilityNeanderthal', 'ProbabilityDenisovan'])
            # add chr prefix
            haps_score['Chromosome'] = 'chr' + haps_score['Chromosome'].astype(str)
            # subset each file by ID so we can run unionbedg for each person separately
            for id in haps_no_score['ID'].unique():
                # get the subset of introgressed haplotypes (no score) for that person
                haps_no_score_id_subset = haps_no_score[haps_no_score['ID']==id]
                # keep specific columns
                haps_no_score_id_subset = haps_no_score_id_subset.loc[:, ['Chromosome', 'Start', 'Stop', 'ID']]
                # save to CSV with ID as the name of the file
                haps_no_score_id_subset.to_csv(f'{output_dir}/tmp/{id}', sep='\t', index=False, header=False)
                # also get the subset of haplotypes WITH SCORE for that person, subsetting specific columns (unionbedg allows 4 cols max)
                haps_score_id_subset = haps_score[haps_score['ID']==id]
                haps_score_id_subset = haps_score_id_subset.loc[:, ['Chromosome', 'Start', 'Stop', 'ProbabilityNeanderthal']]
                # save to CSV
                haps_score_id_subset.to_csv(f'{output_dir}/tmp/{id}_score', sep='\t', index=False, header=False)
                # sort the files
                os.system(f'{BEDTOOLSDIR}/bedtools sort -i {output_dir}/tmp/{id} > {output_dir}/tmp/{id}_sorted')
                os.system(f'{BEDTOOLSDIR}/bedtools sort -i {output_dir}/tmp/{id}_score > {output_dir}/tmp/{id}_score_sorted')
                # merge the files
                os.system(f'{BEDTOOLSDIR}/bedtools merge -i {output_dir}/tmp/{id}_score_sorted -c 4 -o max  > {output_dir}/tmp/{id}_score_merged')
                os.system(f'{BEDTOOLSDIR}/bedtools merge -i {output_dir}/tmp/{id}_sorted -c 4 -o distinct  > {output_dir}/tmp/{id}_merged')
                # run unionbedg to combine scores and haplotypes
                os.system(f'{BEDTOOLSDIR}/bedtools unionbedg -i {output_dir}/tmp/{id}_score_merged {output_dir}/tmp/{id}_merged > {output_dir}/tmp/{id}_final')
                # remove intermmediate files
                os.system(f'rm {output_dir}/tmp/{id} {output_dir}/tmp/*score {output_dir}/tmp/*sorted {output_dir}/tmp/*merged')
            # concatenate each person's data    
            os.system(f'cat {output_dir}/tmp/*_final > {output_dir}/tmp/{pop}')
            # remove intermmediate files
            os.system(f'rm {output_dir}/tmp/*_final')
            # add ancestral group as a separate column
            pop_df = pd.read_csv(f'{output_dir}/tmp/{pop}', sep='\t', low_memory=False, names=['Chromosome', 'Start', 'End', 'Score', 'ID'])
            pop_df['Ancestry'] = pop
            # save pop-specific dfs
            dfs.append(pop_df)

        # then concatenate all the dfs
        final_df = pd.concat(dfs)
        final_df = final_df.sort_values(["Chromosome", "Start", "End"])

        final_df['Chromosome'] = final_df['Chromosome'].str.replace('chr', '')
        
        final_df.to_csv(f'{introgressed_haplotypes_file}', sep='\t', index=False)

        # remove tmp directory
        os.system(f'rm -r {output_dir}/tmp/')

    else:
        
        # read in haplotypes file
        vernot_2016_haps = pd.read_csv(introgressed_haplotypes_file,
                                sep='\t', low_memory=False)
    
        vernot_2016_haps.rename(columns={'Stop': 'End'}, inplace=True)
    
        # make directory to store bedgraphs
        os.system(f'mkdir {output_dir}')
        os.system(f'mkdir {output_dir}/bedgraphs')
        # convert data to bedgraph format
        save_bedgraph(vernot_2016_haps, output_dir+'/bedgraphs')
        # get ancestries as a list and string separated by a space
        ancestries_list = list(vernot_2016_haps['Ancestry'].unique())
        # get list of bedgraph files to be merged in the next step, format is full path separated by " \"
        merged_bg_files = " \\".join([f'{output_dir}'+'/bedgraphs/'+file for file in os.listdir(f'{output_dir}/bedgraphs') if file.startswith("merged_")])
        # get list of ancestries to be used in unionbedg (must be in the same order as the input files)
        ancestries_str = re.sub(f'{output_dir}'+'/bedgraphs/'+r'(merged_|\.bg)', '', merged_bg_files)
        ancestries_str = re.sub(r'(|\.bg)', '', ancestries_str)
        ancestries_str = re.sub(re.escape('\\'), "", ancestries_str)
        # create union file at the population level
        os.system(f'{BEDTOOLSDIR}/bedtools unionbedg -i {merged_bg_files} -header -names {ancestries_str} > {output_dir}/bedgraphs/union_merged.bg')
        # delete intermediate files, keeping only merged files at the ancestral level
        os.system(f'cd {output_dir}/bedgraphs; rm -f sorted*.bg')
        # read in generated union file
        s_final = pd.read_csv(f'{output_dir}/bedgraphs/union_merged.bg', sep='\t')
        # keep only the highest reported score, collapse into one column
        s_final['highest_score'] = s_final[ancestries_list].max(axis=1)
        # drop ancestries
        s_final = s_final.drop(ancestries_list, axis=1)
        # save as a bedgraph
        s_final.to_csv(f'{output_dir}/bedgraphs/union_merged_highest_score.bg', sep='\t', index=False)
    
        # add metadata
        # add IntrogressionSource
        s_final['IntrogressionSource'] = archaic
        # add method name
        s_final['Method']='S*, Vernot 2016'
        # add GenomesSource
        s_final['GenomesSource']=genomes_source
        # add ScoreName
        s_final['ScoreName']='SNPsOnHaplotype'
    
        # save final file
        s_final.to_csv(f'{output_dir}/{output_file_name}', sep='\t', index=False)
        
        # completed
        print(f'Output .bed file is located at {output_dir}/{output_file_name}', sep='')


def sankararaman_2014_raw_to_bed(snp_input_file, hap_input_file, output_dir, output_file_name, genomes_source='1KG and Melanesian samples'):
    """
    This method reads the raw output from browning_2018 (sprime) and cleans it to .bed format with overlapping regions collapsed into one entry at the superpopulation-levels.
    Args:
        snp_input_file is the output of sankararaman_2014_raw_to_bed.sh, containing SNP and haplotype-level data from extended LD and median extended LD files.
        hap_input_file: output of sankararaman_2014_raw_to_bed.sh, generated from haplotype information, has columns Chromosome	Start	Stop	HaplotypeScore	Population
        output_dir: output directory to place the cleaned file 
        output_file_name: output file name to place in the output_dir
        genomes_source: the source that the genomes are from (usually 1000 Genomes Project (1KG) or Simons Genome Diversity Project (SGDP))
    """
    # read in SNP-level file
    sankararaman_2014_hap = pd.read_csv(f'{hap_input_file}', sep='\t', low_memory=False)
    # rename columns
    sankararaman_2014_hap.rename(columns={'Stop': 'End'}, inplace=True)
    # map ancestry
    sankararaman_2014_hap['Ancestry'] = sankararaman_2014_hap['Population'].apply(map_status)
    sankararaman_2014_hap = sankararaman_2014_hap.drop(columns=['Population'])

    # make directory to store bedgraphs
    os.system(f'mkdir {output_dir}')
    os.system(f'mkdir {output_dir}/bedgraphs')
    # convert data to bedgraph format
    save_bedgraph(sankararaman_2014_hap, output_dir+'/bedgraphs')
    # get ancestries as a list and string separated by a space
    ancestries_list = list(sankararaman_2014_hap['Ancestry'].unique())
    # get list of bedgraph files to be merged in the next step, format is full path separated by " \"
    merged_bg_files = " \\".join([f'{output_dir}'+'/bedgraphs/'+file for file in os.listdir(f'{output_dir}/bedgraphs') if file.startswith("merged_")])
    # get list of ancestries to be used in unionbedg (must be in the same order as the input files)
    ancestries_str = re.sub(f'{output_dir}'+'/bedgraphs/'+r'(merged_|\.bg)', '', merged_bg_files)
    ancestries_str = re.sub(r'(|\.bg)', '', ancestries_str)
    ancestries_str = re.sub(re.escape('\\'), "", ancestries_str)
    # create union file at the population level
    os.system(f'{BEDTOOLSDIR}/bedtools unionbedg -i {merged_bg_files} -header -names {ancestries_str} > {output_dir}/bedgraphs/union_merged.bg')
    # delete intermediate files, keeping only merged files at the ancestral level
    os.system(f'cd {output_dir}/bedgraphs; rm -f sorted*.bg')
    # read in generated union file
    crf_final = pd.read_csv(f'{output_dir}/bedgraphs/union_merged.bg', sep='\t')
    # keep only the highest reported score, collapse into one column
    crf_final['highest_score'] = crf_final[ancestries_list].max(axis=1)
    # drop ancestries
    crf_final = crf_final.drop(ancestries_list, axis=1)
    # save as a bedgraph
    crf_final.to_csv(f'{output_dir}/bedgraphs/union_merged_highest_score.bg', sep='\t', index=False)
    # add metadata
    # add IntrogressionSource
    crf_final['IntrogressionSource'] = 'Neanderthal'
    # add method name
    crf_final['Method']='CRF, Sankararaman 2014'
    # add GenomesSource
    crf_final['GenomesSource']=genomes_source
    # add ScoreName
    crf_final['ScoreName']='Neanderthal probability'

    # save final file
    crf_final.to_csv(f'{output_dir}/{output_file_name}', sep='\t', index=False)
    
    # completed
    print(f'Output .bed file is located at {output_dir}/{output_file_name}', sep='')

def sankararaman_2016_raw_to_bed(hap_input_file, output_dir, output_file_name, genomes_source='1KG and Melanesian samples'):
    """
    This method reads the raw output from sankararaman_2016 (CRF) and cleans it to .bed format with overlapping regions collapsed into one entry at the superpopulation-levels.
    Args:
        hap_input_file: output of sankararaman_2016_raw_to_bed.sh, generated from haplotype information
        output_dir: output directory to place the cleaned file 
        output_file_name: output file name to place in the output_dir
        genomes_source: the source that the genomes are from (usually 1000 Genomes Project (1KG) or Simons Genome Diversity Project (SGDP))
    """
    # read in hap-level file
    sankararaman_2016_hap = pd.read_csv(f'{hap_input_file}', sep='\t', low_memory=False)
    # rename columns
    sankararaman_2016_hap.rename(columns={'Stop': 'End'}, inplace=True)
    # map ancestry
    sankararaman_2016_hap['Ancestry'] = sankararaman_2016_hap['Population'].apply(map_status)
    sankararaman_2016_hap = sankararaman_2016_hap.drop(columns=['Population'])

    # make directory to store bedgraphs
    os.system(f'mkdir {output_dir}')
    os.system(f'mkdir {output_dir}/bedgraphs')
    # convert data to bedgraph format
    save_bedgraph(sankararaman_2016_hap, output_dir+'/bedgraphs')
    # get ancestries as a list and string separated by a space
    ancestries_list = list(sankararaman_2016_hap['Ancestry'].unique())
    # get list of bedgraph files to be merged in the next step, format is full path separated by " \"
    merged_bg_files = " \\".join([f'{output_dir}'+'/bedgraphs/'+file for file in os.listdir(f'{output_dir}/bedgraphs') if file.startswith("merged_")])
    # get list of ancestries to be used in unionbedg (must be in the same order as the input files)
    ancestries_str = re.sub(f'{output_dir}'+'/bedgraphs/'+r'(merged_|\.bg)', '', merged_bg_files)
    ancestries_str = re.sub(r'(|\.bg)', '', ancestries_str)
    ancestries_str = re.sub(re.escape('\\'), "", ancestries_str)
    # create union file at the population level
    os.system(f'{BEDTOOLSDIR}/bedtools unionbedg -i {merged_bg_files} -header -names {ancestries_str} > {output_dir}/bedgraphs/union_merged.bg')
    # delete intermediate files, keeping only merged files at the ancestral level
    os.system(f'cd {output_dir}/bedgraphs; rm -f sorted*.bg')
    # read in generated union file
    crf_final = pd.read_csv(f'{output_dir}/bedgraphs/union_merged.bg', sep='\t')
    # keep only the highest reported score, collapse into one column
    crf_final['highest_score'] = crf_final[ancestries_list].max(axis=1)
    # drop ancestries
    crf_final = crf_final.drop(ancestries_list, axis=1)
    # save as a bedgraph
    crf_final.to_csv(f'{output_dir}/bedgraphs/union_merged_highest_score.bg', sep='\t', index=False)
    # add metadata
    # add IntrogressionSource
    crf_final['IntrogressionSource'] = 'Neanderthal'
    # add method name
    crf_final['Method']='CRF, Sankararaman 2016'
    # add GenomesSource
    crf_final['GenomesSource']=genomes_source
    # add ScoreName
    crf_final['ScoreName']='Neanderthal probability'

    # save final file
    crf_final.to_csv(f'{output_dir}/{output_file_name}', sep='\t', index=False)
    
    # completed
    print(f'Output .bed file is located at {output_dir}/{output_file_name}', sep='')


# other helpful functions
def bed_overall_proportion(input_bed_path, background_bed_path, output_coverage_path):
    """
    This method takes a bed file and calculates the proportion of coverage in a genome file
    Args:
        input_bed: input to calculate proportion
        background_bed: background file 
        output_coverage_path: output file path
    """
    # open the overall introgression .bed file
    input_bed = pd.read_csv(input_bed_path, sep='\t', low_memory=False, names=["chrom", "start", "end"])
    # open the background bed file
    background_bed = pd.read_csv(background_bed_path, sep='\t', low_memory=False, names=["chrom", "start", "end"])
    # filter background bed file by the chromosomes included in input .bed file
    filtered_background_bed = background_bed[background_bed['chrom'].isin(input_bed['chrom'])]
    filtered_background_bed.to_csv((os.path.dirname(background_bed_path) + '/filtered_background.bed'), header=False, sep='\t', index=False)
    # run bedtools coverage
    background_bed_dir = os.path.dirname(background_bed_path)
    os.system(f'{BEDTOOLSDIR}/bedtools coverage -a {background_bed_dir}/filtered_background.bed -b {input_bed_path} > {output_coverage_path}')
    # Read the filtered file and compute the sums
    coverage = pd.read_csv(f'{output_coverage_path}', sep='\t', low_memory=False, names=["chrom", "start", "end", "num_background_overlapping_input", "num_input_overlapping_background", "length_background", "fraction_input_overlap_background", 'background'])

    # Compute the ratio
    overall_coverage = np.sum(coverage.num_input_overlapping_background) / np.sum(coverage.length_background)
    
    return(overall_coverage)


def combine_pop_files(list_of_bg_files, output_dir, output_file):
    """
    This function takes list of bg files and concatenates them
    Args:
        list_of_bg_files: list of files to combine
        output_dir: output file path
        output_file: output file name
    """
    pop_concat_df = pd.DataFrame()

    # concatenate files in list
    for file in list_of_bg_files:
        file=pd.read_csv(file, sep='\t', low_memory=False, names=["chrom", "start", "end", 'score'])
        pop_concat_df = pd.concat([pop_concat_df, file])
        pop_concat_df.to_csv(f'{output_dir}/{output_file}', sep='\t', index=False, header=False)
    # sort and merge
    os.system(f'{BEDTOOLSDIR}/bedtools sort -i {output_dir}/{output_file} > {output_dir}/sorted_{output_file}')
    os.system(f'{BEDTOOLSDIR}/bedtools merge -i {output_dir}/sorted_{output_file} -c 4 -o max > {output_dir}/{output_file}')
    # remove intermmediate file
    os.system(f'rm {output_dir}/sorted_{output_file}')

def generate_overlap(input_files_string, output_file_path, boolean, autosomes=False, region=False, input_files_names_str=False, length_threshold=False, x_only=False):
    '''
    This function generates a combined overlap file between multiple .bedgraph files, with chrom/start/end and separate columns denoting the bedgraph values for each input file.
    Args:
        input_files_string: string of merged and sorted input bedgraph file paths to generate union file, separated by space for bedtools unionbedg
        output_file_path: string for the output file path
        boolean: true if you want boolean region overlap, false for keeping original bedgraph scores
        autosomes: true if X chromosome should be exclluded
        input_files_names_str: names to use as headers for each input file in the union file, separated by space        
        region: if overlap is to be computed at the region-level with one basepair overlap -> a shared region
        subset_length: if it contains an integer, all the bedgraphs will contain the same length, cumulatively summing to the specified integer, keeping the highest scores in each file.
        length_threshold: number of base pairs to threshold each input bedgraph file by
        x_only: True if only the X chromosome should be included
    '''
    # increase file limit
    os.system('ulimit -n 5000')
    # if there is a subset_length:
        # subset the bedgraph to fit the length_threshold
    if length_threshold != False:
        for file in input_files_string.split():
            subset_bedgraph_by_score(input_file=f'{file}', 
                                     output_file=f'{file}_subset_highest_score', 
                                     length_threshold=length_threshold, autosomes=autosomes, x_only=x_only)
            # add the new suffix to the input_files_string
            input_files_string = input_files_string.replace(f'{file}', f'{file}_subset_highest_score')
            
    # create union file at the population level
    if input_files_names_str!=False:
        os.system(f'ulimit -n 5000; {BEDTOOLSDIR}/bedtools unionbedg -i {input_files_string} -header -names {input_files_names_str} > {output_file_path}')
    else:
        os.system(f'ulimit -n 5000; {BEDTOOLSDIR}/bedtools unionbedg -i {input_files_string} -header > {output_file_path}')
    # generate and read in union file as dataframe
    overlap = pd.read_csv(f'{output_file_path}', sep='\t', low_memory=False)
    # replace 0s as NA
    overlap.replace(0, np.nan, inplace=True)
    overlap = overlap.copy()
    # add length data
    if region:
        overlap['length'] = 1
    else:
        overlap['length'] = overlap.end - overlap.start
    # List of columns to exclude
    exclude_columns = ['chrom', 'start', 'end', 'length']

    if autosomes==True:
        overlap=overlap.loc[overlap['chrom'] != 'X']
    if x_only == True:
        overlap=overlap.loc[overlap['chrom'] == 'X']
        # remove columns where all values are false (no X)
        overlap=overlap.loc[:, ~(overlap == False).all()]
        
    if boolean==True:
        # Convert non-NA values to True and NA values to False, excluding specified columns
        boolean_overlap = overlap
        boolean_overlap[boolean_overlap.columns.difference(exclude_columns)] = boolean_overlap[boolean_overlap.columns.difference(exclude_columns)].apply(lambda x: pd.notna(x))
        return boolean_overlap
        # save final file
        boolean_overlap.to_csv(f'{output_file_path}', sep='\t', index=False)
    else: 
        return overlap
        # save final file
        overlap.to_csv(f'{output_file_path}', sep='\t', index=False)

def jaccard_similarity(df, cat1, cat2, col_sum_across):
    '''
    df: boolean_overlap dataframe for each category
    '''
    # get the sum of the rows where both columns are true
    intersection = np.sum(df[(df[cat1]== True) & (df[cat2] == True)][col_sum_across])
    # get the sum of the rows where one of or both columns are true
    union = np.sum(df[(df[cat1] | df[cat2] == True)][col_sum_across])
    return intersection / union if union != 0 else 0

def normalized_jaccard_similarity(df, cat1, cat2, col_sum_across):
    '''
    calculate normalized (min set size) jaccard similarities
    Args:
        cat1: str, column 1
        cat2: str, column 2
        col_sum_across: column to sum across, often 'length'
    '''
     # check if both files contain 'x' chromosomes, if not remove X entries
    intersection = np.sum(df[(df[cat1] == True) & (df[cat2] == True)][col_sum_across])
    return intersection / min(np.sum(df[(df[cat1] == True)][col_sum_across]), np.sum(df[(df[cat2] == True)][col_sum_across]))

def plot_jaccard(boolean_df, categories_list, value_column, type='combined', leaf_order=False, region_level=False, combined_title=False, save_fig_path=False, color='Reds'):
    '''
    Makes a heatmap plot of jaccard statistics between maps, using a boolean_df
    Args:
        boolean_df: df with columns for each category with True/False values, and a column for the value to compute by
        categories_list: list of different categories to separate by on the jaccard plot
        value_column: string, the column name to compute values with
        type: str, either 'raw', 'normalized', or 'combined' for raw/normalized/jaccard plots
        leaf_order: order by jaccard values, can be true for automatic ordering or you can put in a list of numbers corresponding to order of labels
        region_level: whether shared fragments are region-level, default is base-pair level
        combined_title: title to include, str
        save_fig_path: str, final file name to save figure
        color: color map to use
    '''
    jaccard_similarities = {}

    # if plotting region-level jaccard, set length to 1 so each region has a weight of 1
    if region_level==True:
        boolean_df[value_column]=1

    # compute jaccard similarities
    for cat1 in categories_list:
        for cat2 in categories_list:
            jaccard_similarities[(cat1, cat2)] = jaccard_similarity(boolean_df, cat1, cat2, value_column)
    heatmap_data = pd.DataFrame(index=categories_list, columns=categories_list)
    for (cat1, cat2), similarity in jaccard_similarities.items():
        heatmap_data.loc[cat1, cat2] = similarity

    # plot if raw jaccard similarities specified
    if type == 'raw':
        # Plot heatmap
        plt.figure(figsize=(10, 8))
        sns.heatmap(heatmap_data.astype(float), annot=True, cmap="Reds", vmin=0, vmax=1)
        plt.title('Jaccard Similarities Heatmap')
                    
    else: # if not raw jaccard specified, compute nomalized jaccard
        normalized_jaccard_similarities = {}
        normalized_jaccard_denominators = {}
        for cat1 in categories_list:
            for cat2 in categories_list:
                normalized_jaccard_similarities[(cat1, cat2)] = normalized_jaccard_similarity(boolean_df, cat1, cat2, value_column)
        normalized_heatmap_data = pd.DataFrame(index=categories_list, columns=categories_list)
        for (cat1, cat2), similarity in normalized_jaccard_similarities.items():
            normalized_heatmap_data.loc[cat1, cat2] = similarity

        if type == 'normalized': # plot norm jaccard if specified
            plt.figure(figsize=(8, 6))
            sns.heatmap(normalized_heatmap_data.astype(float), annot=True, cmap="Reds", vmin=0, vmax=1)
            plt.title('Normalized Jaccard Similarities Heatmap')

        else: # if combined, compute combined jaccard heatmap data
            if leaf_order==True:
                distance_matrix = 1 - heatmap_data
                Z = linkage(distance_matrix, method='complete')

                #combined_jaccard_heatmap
                # Reorder the Jaccard distance matrices according to the leaf order
                leaf_order = leaves_list(Z)
                ordered_samples = distance_matrix.columns[leaf_order]
                # Reorder the entire norm and raw matrices according to the leaf order
                normalized_heatmap_data_reordered = normalized_heatmap_data.loc[ordered_samples, ordered_samples]
                heatmap_data_reordered = heatmap_data.loc[ordered_samples, ordered_samples]
                
                # compute combined matrix
                jaccard_heatmap_combined = heatmap_data_reordered.where(np.triu(np.ones(heatmap_data_reordered.shape)).astype(bool))
                normalized_jaccard_heatmap_combined = normalized_heatmap_data_reordered.where(np.tril(np.ones(normalized_heatmap_data_reordered.shape), k=-1).astype(bool))
                jaccard_heatmap_combined.values[np.tril_indices_from(jaccard_heatmap_combined, k=-1)] = normalized_jaccard_heatmap_combined.values[np.tril_indices_from(normalized_jaccard_heatmap_combined, k=-1)]
            
            # if a leaf order is included as a list
            elif isinstance(leaf_order, list):
                distance_matrix = 1 - heatmap_data
                Z = linkage(distance_matrix, method='complete')

                #combined_jaccard_heatmap
                # Reorder the Jaccard distance matrices according to the leaf order
                leaf_order = leaf_order
                ordered_samples = distance_matrix.columns[leaf_order]
                # Reorder the entire norm and raw matrices according to the leaf order
                normalized_heatmap_data_reordered = normalized_heatmap_data.loc[ordered_samples, ordered_samples]
                heatmap_data_reordered = heatmap_data.loc[ordered_samples, ordered_samples]
                
                # compute combined matrix
                jaccard_heatmap_combined = heatmap_data_reordered.where(np.triu(np.ones(heatmap_data_reordered.shape)).astype(bool))
                normalized_jaccard_heatmap_combined = normalized_heatmap_data_reordered.where(np.tril(np.ones(normalized_heatmap_data_reordered.shape), k=-1).astype(bool))
                jaccard_heatmap_combined.values[np.tril_indices_from(jaccard_heatmap_combined, k=-1)] = normalized_jaccard_heatmap_combined.values[np.tril_indices_from(normalized_jaccard_heatmap_combined, k=-1)]

            else:
                jaccard_heatmap_combined = heatmap_data.where(np.triu(np.ones(heatmap_data.shape)).astype(bool))
                normalized_jaccard_heatmap_combined = normalized_heatmap_data.where(np.tril(np.ones(normalized_heatmap_data.shape), k=-1).astype(bool))
                jaccard_heatmap_combined.values[np.tril_indices_from(jaccard_heatmap_combined, k=-1)] = normalized_jaccard_heatmap_combined.values[np.tril_indices_from(normalized_jaccard_heatmap_combined, k=-1)]
                
            # plot combined jaccard
            plt.figure(figsize=(9, 7))
            heatmap = sns.heatmap(jaccard_heatmap_combined.astype(float), annot=True, cmap=color, vmin=0, vmax=1, annot_kws={"size": 13}, square=True)
            plt.title('Raw Jaccard Similarities', loc='right', fontsize=15)
            if combined_title:
                plt.title(combined_title + '\n', wrap=True, fontsize=20)
            plt.title('Raw Jaccard Similarities', loc='right', fontsize=15)
            plt.xticks(rotation=45, ha='right', fontsize=12)
            plt.yticks(fontsize=12)
            plt.subplots_adjust(bottom=0.15)  # Increase bottom margin for x-axis labels
            heatmap.set_xlabel('Normalized Jaccard Similarities', loc='left', fontsize=15)  # Add x-axis title

    if save_fig_path!=False:
        print('save')
        plt.savefig(f'{save_fig_path}')
        plt.show()

    return heatmap_data_reordered

def plot_cumulative_increase(boolean_overlap_df, fig_save_path=False, region_level=False, title=False, xtitle=False, ytitle=False, rename_methods=False, color_by_pops=False):
    '''
    Plots the cumulative increased in introgressed segments for a boolean_overlap_df
    Args:
        boolean_overlap_df: df with columns for each category with True/False values, and a column for the value to compute by
        color_by_pops: True if you want to include a legend for populations in the 'Ancestry' column
        xtitle: str, x axis title
        ytitle: str, y axis title
        rename_methods: dictionary to rename columns
        region_level: whether shared fragments are region-level, default is base-pair level
        title: title to include, str
        fig_save_path: str, final file name to save figure
    '''

    # drop bed positions to only keep fragment sizes
    boolean_overlap_df = boolean_overlap_df.drop(['chrom', 'start', 'end'], axis=1)
    if rename_methods:
        boolean_overlap_df = boolean_overlap_df.rename(columns=rename_dict)

    # if plotting region-level jaccard, set length to 1 so each region has a weight of 1
    if region_level==True:
        boolean_overlap_df['length']=1
    
    # bar plot of unique values
    # list to store sum of lengths
    length_sums = []
    # list of columns except 'chrom, 'start', 'end', and 'length'
    columns = [col for col in boolean_overlap_df.columns if col not in ['chrom', 'start', 'end', 'length']]
    
    # initialize false across a filter/mask
    cumulative_filter = pd.Series([False] * len(boolean_overlap_df))

    # iterate through each column and add to filter
    for col in columns:

        # get pop from original df
        ancestry = df[df['ID'] == col]['Ancestry']
        # update filter by adding current column
        new_filter = (boolean_overlap_df[col] == True)
        # combine with cumulative filter (keep track of new True values added)
        combined_filter = cumulative_filter | new_filter
        # compute sum of 'length' for newly added True values
        length_sum = boolean_overlap_df[combined_filter]['length'].sum()
        # append result with column name and new length sum
        length_sums.append((col, length_sum, ancestry))
        # update cumulative filter to include current column
        cumulative_filter = combined_filter

    # Sort sums in descending order
    length_sums.sort(key=lambda x: x[1], reverse=False)
    length_sums_df=pd.DataFrame(length_sums)
    length_sums_df.columns=['algorithm', 'count', 'ancestry']
    
    # plot
    plt.figure(figsize=(10, 4))
    
    plt.scatter(length_sums_df.algorithm, length_sums_df['count'], marker='o', color='b')
    
    if xtitle:
        plt.xlabel(xtitle)
    else:
        plt.xlabel('Algorithms (sorted by unique regions added)')
        
    if ytitle:
        plt.xlabel(ytitle)
    else:
        plt.ylabel('Cumulative Unique Genomic Regions')

    if title:
        plt.title(title)
    else:
        plt.title('Cumulative Increase in Unique Genomic Base Pairs by Algorithm')
    plt.xticks([])
    plt.tight_layout(pad=1)

    if fig_save_path==True:
        plt.savefig(fig_save_path)  

def subset_bedgraph_by_score(input_file, length_threshold, output_file=False, autosomes=False, exact=False, x_only=False):
    '''
    for a bedgraph, subset the file into the top regions (totalling a certain length_threshold)
    Args:
        input_file: Pandas DataFrame, containing bedgraph file ('chrom', 'start', 'end', 'score')
        length_threshold: number of bp to threshold
        output_file: file output string
        autosomes=True, auto subset to autosomes only
        exact: the final subset should be exactly the specified length_threshold, if true.
        x_only: True, if x chr only
    '''
    # if df, don't load file
    if isinstance(input_file, pd.DataFrame):
        df=input_file
        df.columns=['chrom', 'start', 'end', 'score']
    else:   
        df = pd.read_csv(input_file, sep='\t', header=None, names=['chrom', 'start', 'end', 'score'])

    if autosomes==True:
        df=df.loc[df['chrom'] != 'X']
    elif x_only==True:
        df=df.loc[df['chrom'] == 'X']

    # calculate length of each region
    df['length'] = df['end'] - df['start']

    # sort DataFrame by score in descending order
    df['score'] = df['score'].astype('float64')
    df_sorted = df.sort_values(by='score', ascending=False).reset_index(drop=True)

    # calculate cumulative sum of lengths
    df_sorted['cumsum_length'] = df_sorted['length'].cumsum()

    # select regions that fully fit in the threshold
    subset = df_sorted[df_sorted['cumsum_length'] <= length_threshold].copy()
    # if first entry > threshold, handle separately
    if len(subset)==0:
        subset = df_sorted.head(1)

    if exact == True:

        # calculate remaining length needed to reach threshold
        remaining_length = length_threshold - subset['cumsum_length'].max()
        
        # if remaining length, keep adding rows threshold is reach
        if remaining_length > 0:
            index = len(subset)  # Start from the next row
            while index < len(df_sorted) and remaining_length > 0:
                next_row = df_sorted.iloc[index].copy()
    
                # if next row is larger than needed, trim
                if next_row['length'] > remaining_length:
                    next_row['end'] = next_row['start'] + remaining_length
                    next_row['length'] = remaining_length
                    remaining_length = 0  # Reached the threshold exactly
                else:
                    # otherwise, add whole row and update remaining length
                    remaining_length -= next_row['length']
    
                # add selected row to subset
                subset = pd.concat([subset, next_row.to_frame().T], ignore_index=True)
                index += 1  # move to next row
        # EXCEPTION: if less than 0 (need to subtract from first entry)
        elif remaining_length < 0:
            # remaining length is neg length_threshold
            remaining_length = -length_threshold
            index = len(subset)  # start from the next row
            # while we don't have an exact match of the length of the first entry/trimmed row
            while remaining_length != 0:
                next_row = df_sorted.iloc[index - 1].copy()
                # if next row is larger than needed, trim it
                if next_row['length'] > -remaining_length:
                    next_row['end'] = next_row['start'] + -remaining_length
                    next_row['length'] = -remaining_length
                    remaining_length = 0
    
                # add selected row to the subset
                subset = next_row.to_frame().T
                index += 1  # Move to the next row
    
        if subset['length'].sum()!=length_threshold:
            print('Error, total length does not equal the threshold:', subset['length'].sum())

    # drop helper columns
    subset = subset.drop(columns=['length', 'cumsum_length'])

    # sort the file
    subset = subset.sort_values(["chrom", "start", "end"])
    # convert to integers if autosomes only
    if 'X' not in subset['chrom'].values:
        subset['chrom']=subset['chrom'].astype('int')
    else:
        # remove trailing zeros
        subset['chrom']=subset['chrom'].astype(str).str.rstrip('0').str.rstrip('.')
    subset['start']=subset['start'].astype('int')
    subset['end']=subset['end'].astype('int')

    # if output is specified
    if output_file != False:
        # aave the subset to a new bedgraph file
        subset.to_csv(output_file, sep='\t', header=False, index=False)
    else:
        return subset


def individual_jaccard_comparison(dfs_dict, temp_dir, output_jaccards_file_path, autosomes_only, x_only, save_shared_individuals_bed=False, length_threshold=False, return_heatmap=False):
    '''
    computed jaccard across shared individuals between methods.

    dfs_dict: dictionary containing algorithm names as keys and the corresponding algorithm's individual-level df
        individual level DFs columns needed:
                                Chromosome	
                                Start	
                                End	
                                Score	
                                Population	
                                Ancestry	
                                ID: person ID, Sample_ID(Aliases) from sgdp metadata and individual IDs for 1000 genomes
                                ScoreName	
                                IntrogressionSource	
                                Method	
                                GenomesSource
    temp_dir: temporary dir to put data in
    output_jaccards_file_path: string, where to put the final .tsv file
    autosomes_only: subset to only look at autosomes, must be used with x_only
    x_only: subset to only look at X chromosome
    save_shared_individuals_bed: save the bed file with shared individuals between methods, should be a string for the file path
    length_threshold: total length to subset a bedgraph by
    return_heatmap: True if you want to return the heatmap
    '''
    # report number of IDs shared between methods
    # get list of different algorithms in comparison
    num_algorithms = len(dfs_dict)
    algorithms = list(dfs_dict.keys())

    # make combined dataframe
    individual_concat = pd.concat([df.assign(Algorithm=k) for k,df in dfs_dict.items()])
    # keep autosomes only
    if autosomes_only == True:
        individual_concat = individual_concat[individual_concat['Chromosome'] != 'X']
    if x_only == True:
        individual_concat = individual_concat[individual_concat['Chromosome'] == 'X']
        
    # group by ID and count unique values in 'Value'
    id_value_counts = individual_concat.groupby('ID')['Algorithm'].nunique()
    # filter IDs in all individual level dataframes
    shared_individuals = id_value_counts[id_value_counts == num_algorithms].index.tolist()    
    print(algorithms, ":", len(shared_individuals), 'individuals shared, ', shared_individuals)

    shared_individuals_bed = individual_concat[individual_concat['ID'].isin(shared_individuals)]
    if save_shared_individuals_bed != False:
        shared_individuals_bed.to_csv(f'{save_shared_individuals_bed}', sep='\t', index=False)

    # list to store final heatmap data
    final_jaccards = []
    final_normalized_jaccards = []
    # list to store person IDs
    persons = []
    person_ids = [] # this one also stores IDs, in the order that they are in the final tsv file

    # for each ID, calculate Jaccard
    for person in shared_individuals:
        
        # store jaccards in dict
        jaccard_similarities = {}
        normalized_jaccard_similarities = {}

        # if there is a length_threshold:
        if length_threshold != False:
            # compute the min length between individuals to use
            # compute length
            shared_individuals_bed['length'] = shared_individuals_bed['End'] - shared_individuals_bed['Start']
            # empty dict to hold proportions
            proportions = {}
            # for each algorithm, subset df by the method
            for method in shared_individuals_bed.Algorithm.unique():
                method_subset = shared_individuals_bed[shared_individuals_bed['Algorithm'] == method]
                # compute the minumum amount of introgression across people
                min_sum = method_subset.groupby('ID')['length'].sum().min()
                # add the minimum amount to the dict
                proportions[method] = min_sum
                # get the smallest value to use as length_threshold
                length_threshold = min(proportions.values())
    
        # separate each person into their own dataframe as bedgraphs, storing it in a dictionary with keys as the alg name and value as the individual's bedgraph
        individual_bedgraphs_dict={}
        for algorithm_name, df in dfs_dict.items():
            individual_bedgraph=df[df['ID']==person][['Chromosome', 'Start', 'End', 'Score']]

            # keep autosomes only
            if autosomes_only == True:
                individual_bedgraph = individual_bedgraph[individual_bedgraph['Chromosome'] != 'X']
            if x_only == True:
                individual_bedgraph = individual_bedgraph[individual_bedgraph['Chromosome'] == 'X']

            # for ArchIE, handle the semicolon scores
            # keep the highest score reported for a haplotype
            if individual_bedgraph['Score'].astype(str).str.contains(';').any():
                individual_bedgraph['Score'] = individual_bedgraph['Score'].apply(lambda x: max(map(float, x.split(';'))) if isinstance(x, str) else x)
        
            # temporarily save the dfs to run bedtools on them
            individual_bedgraph.to_csv(f'{temp_dir}/{algorithm_name}_{person}', sep='\t', index=False, header=False)

            # sort bedgraphs
            os.system(f'{BEDTOOLSDIR}/bedtools sort -i {temp_dir}/{algorithm_name}_{person} > {temp_dir}/sorted_{algorithm_name}_{person}')
    
            # merge bedgraphs
            os.system(f'{BEDTOOLSDIR}/bedtools merge -i {temp_dir}/sorted_{algorithm_name}_{person} -c 4 -o max > {temp_dir}/{algorithm_name}_{person}')

            # if there is a length_threshold:
            # subset the bedgraph to fit the length_threshold
            if length_threshold != False:
                subset_bedgraph_by_score(input_file=f'{temp_dir}/{algorithm_name}_{person}', 
                                         output_file=f'{temp_dir}/{algorithm_name}_{person}', 
                                         length_threshold=length_threshold)
        
        # get string for files to include
        files_list = []
        for algorithm_name in algorithms:
            files_list.append(f'{temp_dir}/{algorithm_name}_{person}')
        # flatten the list and add " \" as the separator for bash
        files_string = ' '.join(files_list)
        algorithms_string = ' '.join(algorithms)
        
        # create overlap file
        overlap=generate_overlap(input_files_string=f'{files_string}',
                           input_files_names_str=f'{algorithms_string}',
                           output_file_path=f'{temp_dir}/{person}',
                           boolean=True)

        # compute jaccard stats
        for cat1 in algorithms:
            for cat2 in algorithms:
                jaccard_similarities[(cat1, cat2)] = jaccard_similarity(overlap, cat1, cat2, 'length')
        # create heatmap df
        heatmap_data = pd.DataFrame(index=algorithms, columns=algorithms)
        # save jaccards to heatmap df
        for (cat1, cat2), similarity in jaccard_similarities.items():
            heatmap_data.loc[cat1, cat2] = similarity
    
        final_jaccards.append(heatmap_data)
        
        # add ID to running list to add to final dataframe
        person_ids.append(person)
        
        #remove temp files
        os.system(f'rm {temp_dir}/{person}')
        
        for algorithm_name in algorithms:
            os.system(f'rm {temp_dir}/{algorithm_name}_{person} \
                           {temp_dir}/sorted_{algorithm_name}_{person}')

        persons.append(person)

    # flatten every heatmap dataframe -> single column
    flattened_columns = []
    for i, df in enumerate(final_jaccards):
        # melt dataframe to long format
        melted_df = df.reset_index().melt(id_vars='index', var_name='column', value_name='Jaccard')
        # create combined identifier for the pairs (method1:method2)
        melted_df['Pair'] = melted_df['index'] + ':' + melted_df['column']
        # set 'pair' column as index
        melted_df.set_index('Pair', inplace=True)
        # add corresponding person ID
        melted_df['ID'] = persons[i]  
        # append melted DataFrame
        flattened_columns.append(melted_df[['ID', 'Jaccard']])
    
    # combine all into final dataframe
    flat_df = pd.concat(flattened_columns).reset_index().rename(columns={'index': 'Pair'})
    # save final dataframe
    flat_df.to_csv(f'{output_jaccards_file_path}', sep='\t', index=False)

    if return_heatmap==True:
        return final_jaccards
    else:
        return flat_df

def individual_jaccard_comparison_heatmap(tsv, title, color, save_plot_path, figsize_tuple=(8, 6)):
    '''
    plots the lower triangle heatmap of the comparison of individual level jaccards between methods
    Args:
        df: str, tsv path for a dataframe containing columns 'Pair' and 'Jaccard', where each row contains pair, a string containing 2 method names separated by a '-', and the corresponding jaccard
        title: str, title for the jaccard plot
        color: str, color map parameter for the plot, cmaps in matplot lib https://matplotlib.org/stable/users/explain/colors/colormaps.html
        save_plot_path: str, file path to save the final plot
        figsize_tuple: tuple, adjust figure size
    '''
    
    if isinstance(tsv, pd.DataFrame):
        df=tsv
    else:
        df=pd.read_csv(f'{tsv}', sep='\t', low_memory=False)
    
    # sort pairs alphabetically
    df['Sorted_Pair'] = df['Pair'].apply(lambda x: ':'.join(sorted(x.split(':'))))

    # compute average between pairs
    df_avg = df.groupby('Sorted_Pair', as_index=False).agg({'Jaccard': 'mean'})
    
    # create a pivot table for the heatmap, splitting the pair into two columns
    df_avg['Name1'], df_avg['Name2'] = zip(*df_avg['Sorted_Pair'].map(lambda x: x.split(':')))
    
    # create pivot table with Jaccard values
    heatmap_data = df_avg.pivot(index='Name1', columns='Name2', values='Jaccard')
    
    # fill the other half to make a symmetric matrix
    heatmap_data = heatmap_data.combine_first(heatmap_data.T)

    # remove diagonal
    np.fill_diagonal(heatmap_data.values, np.nan)
    
    # create a mask for the upper triangle and diagonal
    mask = np.zeros_like(heatmap_data, dtype=bool)
    np.fill_diagonal(mask, True)  # mask diagonal
    mask = np.triu(np.ones_like(mask, dtype=bool))
    
    # figure size
    plt.figure(figsize=figsize_tuple)
    # plot heatmap with mask to hide upper triangle
    ax=sns.heatmap(heatmap_data, annot=True, cmap=f'{color}', linewidths=.5, mask=mask, vmin=0, vmax=1)

    # get current labels
    x_labels = ax.get_xticklabels()
    y_labels = ax.get_yticklabels() 

    # set new labels, excluding first y-label and last x-label
    new_x_labels = [label.get_text() for label in x_labels[:-1]]
    new_y_labels = [label.get_text() for label in y_labels[1:]]
    
    # update ticks and labels on the heatmap
    ax.set_xticks(ax.get_xticks()[:-1])
    ax.set_yticks(ax.get_yticks()[1:])
    
    # apply new labels
    ax.set_xticklabels(new_x_labels)
    ax.set_yticklabels(new_y_labels)

    # add title
    plt.title(f'{title}')
    # rotate x
    plt.xticks(rotation=45, ha='right')

    # expand final plot
    plt.tight_layout()
    
    plt.savefig(f'{save_plot_path}')

    plt.show()


def individual_jaccard_comparison_dot_heatmap(heatmap_data, tsv, color, save_plot_path, figsize_tuple=(8, 6), text_position=0.2, sd_sizes=False, title=False, legend_title=False):
    '''
    Plots a lower triangle dot heatmap showing the mean and SD of individual-level Jaccard comparisons between methods.
    
    Parameters:
    tsv (str or pd.DataFrame): Path to TSV file or dataframe containing columns 'Pair' and 'Jaccard'.
    title (str): Title for the plot.
    color (str): Colormap for SD representation.
    save_plot_path (str): Path to save the plot.
    figsize_tuple (tuple): Figure size.
    text_position: distance from the top to plot text
    sd_sizes: list of values for SD sizes in the legend
    legend_title (str): optional title for the legend
    '''
    
    if len(heatmap_data) != 1:
        # compute the average Jaccards between the list of jaccard dataframes (in heatmap_data) if there is more than one heatmap (person)
        df = pd.concat(heatmap_data).groupby(level=0).mean()
    else:
        df = heatmap_data[0]
    
    # reorder and match up columns
    df_reordered = df.reindex(columns=df.columns, index=df.columns)
    # compute distance and linkage matrices
    distance_matrix = 1 - df_reordered
    # linkage matrices
    Z = linkage(distance_matrix, method='complete')
    
    # reorder jaccards according to leaf order
    leaf_order = leaves_list(Z)
    ordered_samples = distance_matrix.columns[leaf_order]
    
    # load data
    if isinstance(tsv, pd.DataFrame):
        df = tsv
    else:
        df = pd.read_csv(tsv, sep='\t', low_memory=False)
    
    # sort pairs alphabetically
    df['Sorted_Pair'] = df['Pair'].apply(lambda x: ':'.join(sorted(x.split(':'))))
    
    # compute mean and SD
    df_stats = df.groupby('Sorted_Pair', as_index=False).agg({'Jaccard': ['mean', 'std']})
    df_stats.columns = ['Sorted_Pair', 'Mean_Jaccard', 'SD_Jaccard']
    
    # split into two columns
    df_stats['Name1'], df_stats['Name2'] = zip(*df_stats['Sorted_Pair'].map(lambda x: x.split(':')))
    
    # pivot for heatmap,reorder by distance
    jaccard_data = df_stats.pivot(index='Name1', columns='Name2', values='Mean_Jaccard').loc[ordered_samples, ordered_samples]
    sd_data = df_stats.pivot(index='Name1', columns='Name2', values='SD_Jaccard').loc[ordered_samples, ordered_samples]
    
    # fill NaNs with symmetric values
    jaccard_data = jaccard_data.where(~jaccard_data.isna(), jaccard_data.T)
    sd_data = sd_data.where(~sd_data.isna(), sd_data.T)
    # mask upper triangle
    jaccard_data = jaccard_data.mask(np.triu(np.ones(jaccard_data.shape), k=1).astype(bool))
    sd_data = sd_data.mask(np.triu(np.ones(sd_data.shape), k=1).astype(bool))
    
    # create mask for upper triangle
    mask = np.triu(np.ones_like(jaccard_data, dtype=bool))
    
    # create figure
    fig, ax = plt.subplots(figsize=(8, 6))
    
    # plot heatmap
    if legend_title:
        heatmap = sns.heatmap(jaccard_data, cmap="Blues", linewidths=0.5, cbar_kws={'label': legend_title}, ax=ax, annot=False, vmin=0, vmax=1)
    else:
        heatmap = sns.heatmap(jaccard_data, cmap="Blues", linewidths=0.5, cbar_kws={'label': 'Mean Jaccard'}, ax=ax, annot=False, vmin=0, vmax=1)
    
    # overlay mean values
    for i in range(jaccard_data.shape[0]):
        for j in range(jaccard_data.shape[1]):
            if not np.isnan(jaccard_data.iloc[i, j]):  # Ensure valid values
                ax.text(j + 0.5, i + text_position, f"{jaccard_data.iloc[i, j]:.2f}", ha='center', va='center', color='black', fontsize=12)
    
    # overlay scatter plot for SD Jaccard (move downward)
    for i in range(sd_data.shape[0]):
        for j in range(sd_data.shape[1]):
            if not np.isnan(sd_data.iloc[i, j]):  # Ensure valid SD values
                ax.scatter(j + 0.5, i + 0.6, s=sd_data.iloc[i, j] * 5000, color='black', alpha=0.6, edgecolors='white')
    
    # create dot size legend if specified
    if sd_sizes:
        sizes = sd_sizes
    else:
        sizes = [sd_data[sd_data != 0].min().min(), sd_data[sd_data != 0].max().max()]  # SD values for legend
    legend_handles = [plt.scatter([], [], s=s * 5000, color='black', alpha=0.6, edgecolors='white', label=f'SD: {s:.2f}') for s in sizes]
    
    # add legends
    legend1 = plt.legend(handles=legend_handles, title="Jaccard SD", loc="upper right", frameon=True, fontsize=10)
    ax.add_artist(legend1)  # Add dot legend separately
    
    # labels and title
    if title:
        plt.title(title)
    # rotate x
    plt.xticks(rotation=45, ha='right')
    
    # aave and show plot
    plt.tight_layout()
    plt.savefig(save_plot_path)
    plt.show()