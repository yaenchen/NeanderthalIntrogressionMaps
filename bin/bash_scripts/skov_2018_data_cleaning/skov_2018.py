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


tools.skov_2018_raw_to_bed(input_folder=f'{HOMEDIR}/introgression_methods/data/skov_2018', 
                           output_dir=f'{HOMEDIR}/introgression_methods/cleaned/introgression_tools/skov_2018/', 
                           liftover_path='/wynton/group/capra/bin/liftOver',
                           output_file_name='nean_skov2018_frag.bed', 
                           archaic='Neanderthal',
                           individual_level=True,
                           scratch_dir=sys.argv[1]) # scratch from command line