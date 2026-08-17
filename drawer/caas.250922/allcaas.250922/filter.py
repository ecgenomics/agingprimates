import pandas as pd
import sys

with open(sys.argv[1], "r") as handle:
    df = pd.read_csv(handle, sep="\t", header=0)

filtered = df[
        (df["Pvalue"] < 0.05)]

    # Save to new file (append .filtered.tsv)
outname = sys.argv[1] + ".filtered.tsv"
filtered.to_csv(outname, sep="\t", index=False)

print(f"{sys.argv[1]}: {len(filtered)} rows kept → saved in {outname}")