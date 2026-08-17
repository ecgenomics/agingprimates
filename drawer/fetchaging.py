#!/usr/bin/env python

"""
AGING IN PRIMATES: simplified data processing script

This script computes the Longevity Quotient (LQ) for primate species
based on their maximum lifespan and body mass.

Input:
- Trait definitions (MLS and BM)
- Phenomic dataset with trait values

Output:
- Tab-separated file with LQ values and basic traits per species
"""

import pandas as pd
import os

### INPUT FILES

bm_traits = "inputs/bm.traits.manual.tab"
mls_traits = "inputs/mls.traits.manual.tab"
phenomic_data = "inputs/nhp.phenomic.dataset.csv"
genomic_data = "inputs/tmb.sample.count.tab"

### OUTPUT FILE

lqdf_out = "dataset/lqdf.tab"

### FUNCTION: compute_LQ

def compute_LQ(lifespan_years, body_mass_g):
    """
    Compute Longevity Quotient (LQ) based on body mass.
    Equation from de Magalhaes et al. (2007):
        LQ = MLS / (6.47 * M^0.189)
    """
    try:
        a = 6.47
        b = 0.189
        expected_lifespan = a * (body_mass_g ** b)
        return lifespan_years / expected_lifespan
    except Exception as e:
        print(f"⚠️ Error in LQ computation: {e}")
        return None


### FUNCTION: load_traits

def load_traits(filein):
    class traitinfo():
        def __init__(self):
            self.codes = []
            self.names = []
            self.d = {}

    z = traitinfo()
    with open(filein, "r") as h:
        l = h.read().splitlines()
    z.codes = [x.split("\t")[0] for x in l]
    z.names = [x.split("\t")[1] for x in l]
    z.d = dict(zip(z.codes, z.names))
    return z

### PREPARE OUTPUT FOLDER

os.makedirs("dataset", exist_ok=True)

### LOAD GENETIC INFO

with open(genomic_data) as h:
    l = h.read().splitlines()
    genome_availability = {x.split("\t")[0] : int(x.split("\t")[1]) for x in l}

### LOAD TRAIT INFO

bm = load_traits(bm_traits)
mls = load_traits(mls_traits)

### LOAD PHENOMIC DATA

phenomic_df = pd.read_csv(phenomic_data)

### SELECT AND RENAME RELEVANT COLUMNS

selected_columns = ['GroupName'] + mls.codes + bm.codes
lqdf = phenomic_df[selected_columns]

# Rename columns to human-readable names
lqdf = lqdf.rename(columns={**mls.d, **bm.d}, errors='ignore')

# Parse species name
split_cols = lqdf['GroupName'].str.split('_', expand=True)
lqdf['Superfamily'] = split_cols[0]
lqdf['Family'] = split_cols[1]
lqdf['Species'] = split_cols[2] + "_" + split_cols[3]

### COMPUTE LQ

# Assume lifespan is in years and mass is in kg (convert to grams)
lqdf['LQ'] = lqdf.apply(
    lambda row: compute_LQ(row['Maximum.lifespan..y.'], row['BodyMass_kg'] * 1000),
    axis=1
)

lqdf["GenomeSamples"] = lqdf["Species"].map(genome_availability).fillna(0).astype(int)

# Drop rows with missing LQ
lqdf = lqdf.dropna(subset=['LQ'])

# Compute z-score of LQ
lq_mean = lqdf['LQ'].mean()
lq_std = lqdf['LQ'].std()
lqdf['LQ_zscore'] = (lqdf['LQ'] - lq_mean) / lq_std

### BUILD COLUMNS

final_columns = ['Superfamily', 'Family', 'Species', 'Maximum.lifespan..y.', 'BodyMass_kg', 'LQ', 'LQ_zscore', "GenomeSamples"]
lqdf = lqdf[final_columns]
lqdf = lqdf[lqdf["GenomeSamples"] > 0]

### ZSCORE PER FAMILY

# Calcola media e std del LQ per ogni famiglia
family_stats = lqdf.groupby('Family')['LQ'].agg(['mean', 'std']).rename(
    columns={'mean': 'LQ_mean_family', 'std': 'LQ_std_family'}
)

# Unisci queste statistiche al DataFrame originale
lqdf = lqdf.merge(family_stats, left_on='Family', right_index=True)

# Calcola lo z-score rispetto alla media familiare
lqdf['LQ_zscore_family'] = (lqdf['LQ'] - lqdf['LQ_mean_family']) / lqdf['LQ_std_family']


### LOCKING USELESS COLUMNS
final_columns = ['Superfamily', 'Family', 'Species', 'Maximum.lifespan..y.', 'BodyMass_kg', 'LQ', 'LQ_mean_family']
lqdf = lqdf[final_columns]
# Species with exceptional longevity (LQ > 1.5)
lqdf['LongLivedCFG'] = (lqdf['LQ'] > 1.3).astype(int)

# Species with much shorter-than-expected lifespan (LQ < 0.5)
lqdf['ShortLivedCFG'] = (lqdf['LQ'] < 0.7).astype(int)

### EXPORT RESULT

# Export: Species vs LongLivedCFG
lqdf[['Species', 'LongLivedCFG']].to_csv(
    "dataset/species_longlivedcfg.tab", sep="\t", header=False, index=False
)

# Export: Species vs ShortLivedCFG
lqdf[['Species', 'ShortLivedCFG']].to_csv(
    "dataset/species_shortlivedcfg.tab", sep="\t", header=False, index=False
)

lqdf.to_csv(lqdf_out, sep="\t", index=False)