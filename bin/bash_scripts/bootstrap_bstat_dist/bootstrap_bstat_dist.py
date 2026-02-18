#!/bin/env python

import gzip
import pandas as pd
pd.options.mode.chained_assignment = None
import numpy as np
import matplotlib.pyplot as plt
# set font to helvetica for the manuscript
import matplotlib as mpl
from matplotlib import pyplot
# use custom file with helvetica as default
import matplotlib.font_manager as fm
font_path = "/wynton/home/capra/ychen39/.local/lib/python3.11/site-packages/matplotlib/mpl-data/fonts/ttf/Helvetica.ttf"
fm.fontManager.addfont(font_path)
mpl.rcParams['font.family'] = 'Helvetica'
import seaborn as sns
import pyranges as pr
import pyBigWig
import os
import pybedtools
BEDTOOLSDIR = '/wynton/home/cbi/shared/software/CBI/bedtools2-2.31.1/bin'
pybedtools.helpers.set_bedtools_path(path=f'{BEDTOOLSDIR}')
import tarfile
import io
import sklearn
from sklearn.preprocessing import LabelEncoder
import re
import upsetplot
from upsetplot import UpSet, from_indicators, plot
from scipy.cluster.hierarchy import dendrogram, linkage, leaves_list
from collections import defaultdict
import itertools
import joypy
from pathlib import Path
from scipy import stats
from scipy.stats import mannwhitneyu
import pickle

## WYNTON
HOMEDIR='/wynton/home/capra/ychen39/'

# renaming dict
rename_dict = {'argweaverd': 'ARGweaver-D', 'archaicseeker': 'ArchaicSeeker2', 'archie': 'ArchIE',
                                   'sprime': 'Sprime', 'sarge': 'SARGE', 'vernot_2016':'S*', 'sankararaman_2014': 'CRF14',
                                 'sankararaman_2016_1':'CRF16 (1)', 'sankararaman_2016_2':'CRF16 (2)', 'skov_2020':'hmmix20',
                                 'steinruecken_2018':'DICAL-ADMIX', 'vernot_2016_extendedLD': 'S* (extended LD)', 
                                 'vernot_2016_medianextendedLD': 'S* (median extended LD)', 'ibdmix24':'IBDMix', 'iasi_2024':'admixfrog', 'skov_2018':'hmmix18'} # added after reviews

rename_dict_dir = {'hubisz_2020': 'ARGweaver-D', 'yuan_2021': 'ArchaicSeeker2', 'durvasula_2019': 'ArchIE',
                                   'browning_2018_updated_filtering': 'Sprime', 'schaefer_2021': 'SARGE', 'vernot_2016':'S*', 'sankararaman_2014': 'CRF14',
                                 'sankararaman_2016_1':'CRF16 (1)', 'sankararaman_2016_2':'CRF16 (2)', 'skov_2020':'hmmix20',
                                 'steinruecken_2018':'DICAL-ADMIX', 'li_2024':'IBDMix', 'iasi_2024':'admixfrog', 'skov_2018':'hmmix18'} # added after reviews

methods_dirs = ['yuan_2021', 'durvasula_2019', 'steinruecken_2018', 'skov_2020', 'schaefer_2021', 
                'hubisz_2020', 'browning_2018_updated_filtering', 'vernot_2016', 'sankararaman_2014', 'sankararaman_2016_1', 'sankararaman_2016_2',
                'li_2024', 'iasi_2024', 'skov_2018'] # added after reviews

def bootstrap_mean_ci(data, n_bootstrap=1000, ci=95):
    """returns bootstrapped mean and confidence interval bounds."""
    means = []
    for _ in range(n_bootstrap):
        sample = np.random.choice(data, size=len(data), replace=True)
        means.append(np.mean(sample))
    lower = np.percentile(means, (100 - ci) / 2)
    upper = np.percentile(means, 100 - (100 - ci) / 2)
    return np.mean(means), lower, upper


# open precomputed bstat for each window
with open(f'{HOMEDIR}/introgression_methods/data/bstat_by_map_500bp.pkl', 'rb') as f:
    combined = pickle.load(f)


# group, bootstrap
results = []

for (feature, source), group in combined.groupby(["Feature", "Source"]):
    mean, ci_lower, ci_upper = bootstrap_mean_ci(group["Value"].values)
    results.append({
        "Feature": feature,
        "Source": source,
        "Mean": mean,
        "CI_lower": ci_lower,
        "CI_upper": ci_upper
    })

ci_df = pd.DataFrame(results)


ci_df.to_csv(f'{HOMEDIR}/introgression_methods/data/bstat_by_map_500bp_bootstrapped.tsv', sep='\t', index=False, header=False)

# plot averages
plt.figure(figsize=(10, 6))
sns.barplot(
    data=ci_df,
    x="Feature",
    y="Mean",
    hue="Source",
    capsize=0.1,
    palette="Set2"
)

# add error bars manually based on bootstrapping
for i, row in ci_df.iterrows():
    x = list(ci_df["Feature"].unique()).index(row["Feature"])
    offset = -0.2 if row["Source"] == ci_df["Source"].unique()[0] else 0.2
    plt.errorbar(
        x + offset,
        row["Mean"],
        yerr=[[row["Mean"] - row["CI_lower"]], [row["CI_upper"] - row["Mean"]]],
        fmt='none',
        c='black',
        capsize=5
    )

plt.ylabel("Bootstrapped Mean Value")
plt.title("Mean Value by Feature and Source with 95% CI")
plt.tight_layout()
pyplot.savefig(f'{HOMEDIR}/introgression_methods/figures/introgression_tools/bstatscores_bootstrap.svg', bbox_inches='tight')
plt.show()