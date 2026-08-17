#!/usr/bin/env python3

### Aging in primates.
### Binarisation of contrasts. Create the table for decision on contrast binarisation.

import pandas as pd
from datetime import datetime
import os
from itertools import product


# Inputs

pairs = "pairs_overview_other_traits.tsv"    # Structure: Trait\tTop\tBottom
lqtable = "../dataset/lqdf.tab"              # The longevity quotient table

# Step 1: read the inputs and create the dictionaries

#   1.1 Read the tables
pairs_df = pd.read_table(pairs, sep="\t")
lq_df = pd.read_table(lqtable, sep="\t")

#   1.2 Creation of dictionaries
species2family = dict(zip(lq_df['Species'], lq_df['Family']))
species2LQ = dict(zip(lq_df['Species'], lq_df['LQ']))
family2LQ = dict(zip(lq_df['Family'], lq_df['LQ_mean_family']))

# Step 2: Create the table

#   2.1 Add family and reorder
ot = pairs_df[["Top", "Bottom"]]
ot["Family"] = ot["Top"].map(species2family)
ot = ot[["Family"] + [col for col in ot.columns if col != "Family"]]

#   2.2 Add LQ values

#       2.2.1 Family Average (LQ.fam.avg)
ot["LQ.fam.avg"] = ot["Family"].map(family2LQ)

#       2.2.2 Top LQ (LQ.top)
ot["LQ.top"] = ot["Top"].map(species2LQ)

#       2.2.2 Bottom LQ (LQ.bottom)
ot["LQ.bottom"] = ot["Bottom"].map(species2LQ)

#   2.3 Calculate LQ variation to family (LQ.var.fam)
ot["LQ.top.diff.fam"] = abs(ot["LQ.top"] - ot["LQ.fam.avg"])
ot["LQ.bottom.diff.fam"] = abs(ot["LQ.bottom"] - ot["LQ.fam.avg"])

ot["Flag"] = (ot["LQ.top.diff.fam"] > 0.2) & (ot["LQ.bottom.diff.fam"] > 0.2)

#   sort the table
ot = ot.sort_values(by='Family', ascending=True)


# Step 3: Print the table

timestamp = datetime.now().strftime('%y%m%d')

ot.to_csv(f'table.{timestamp}.tsv', sep='\t', index=False)

# Step 4: Build the cfg files.

#   4.1: Folder creation

foldername = "cfg." + timestamp

try:
    os.mkdir(foldername)
except:
    os.system("rm -r " + foldername)
    os.mkdir(foldername)

fam2top = (ot.groupby("Family")["Top"].apply(list).to_dict())
fam2bottom = (ot.groupby("Family")["Bottom"].apply(list).to_dict())

cfg_generic_string = ""

for x in fam2top:
    if len(fam2top[x]) > 1:
        top = "@t"+x
    else:
        top = fam2top[x][0]
    cfg_generic_string = cfg_generic_string + top + "\t1\n"

for x in fam2bottom:
    if len(fam2bottom[x]) > 1:
        bottom = "@b"+x
    else:
        bottom = fam2bottom[x][0]
    cfg_generic_string = cfg_generic_string + bottom + "\t0\n"

# Replace with multiple species.

multifam = [x for x in fam2top if len(fam2top[x]) > 1]

indexes = (list(range(len(fam2top[x]))) for x in multifam)
combs = list(product(*indexes))

counter = 0
for x in combs:
    counter += 1
    out = open(foldername + "/6v6.hyp" + str(counter) +  ".cfg", "w")
    k = cfg_generic_string
    for i, val in enumerate(x):
        fam = multifam[i]
        k = k.replace("@t"+fam, fam2top[fam][val]).replace("@b"+fam, fam2bottom[fam][val])
    out.write(k)
    out.close()