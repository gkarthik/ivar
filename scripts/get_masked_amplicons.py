#!/usr/local/bin/python3.5
import sys, os
import pandas as pd
from scipy.stats import fisher_exact as fe
from statsmodels.sandbox.stats.multicomp import multipletests

prefix = sys.argv[1]
freq = float(sys.argv[2])
df_bed_paths = sys.argv[3:]
df_paths = df_bed_paths[0::2]
bed_paths = df_bed_paths[1::2]

out_paths = []
for i in range(len(df_paths)):
    _ = os.path.basename(df_paths[i])
    _ = _.split(".")[0]
    out_paths.append(_)

beds = []
for i in bed_paths:
    beds.append(pd.read_table(i, sep="\t", names=["Region", "Start", "End", "Name", "Score", "Strand"]))

vdfs = []
for _i,i in enumerate(df_paths):
    _ = pd.read_table(i, names=["POS", "REF", "ALT", "AD", "REV", "DP", "QUAL"], skiprows = 1)
    vdfs.append(_)

# Fisher's Exact Test
#           | AD | DP  |
# Variant   |    |     |
# Threshold | 3  | 100 |
for vdf in vdfs:
    pvals = vdf.apply(lambda x: fe([[x["AD"], x["DP"]], [(freq/100)*x["DP"], x["DP"]]], "greater"), axis = 1)
    vdf["threshold_"+str(freq)+"%_pval"] = multipletests([i[1] for i in pvals], method="fdr_bh")[1]
    vdf["threshold_"+str(freq)+"%_oddsratio"] = [i[0] for i in pvals]

threshold = freq
col = "threshold_"+str(threshold)+"%"
pval_threshold = 0.05
_ = vdfs
# _ = [i[i[col+"_pval"]<=pval_threshold] for i in _]
df = _[0]
for _i, i in enumerate(_[1:]):
    df = df.merge(i, how='inner', on=['POS', 'REF', 'ALT'], suffixes = ("_0", "_"+str(_i+1)))

cols = df.columns[df.columns.str.match(r"\b"+col+"\b*_pval")]
df = df.ix[df[cols].apply(lambda x: any([i<= pval_threshold for i in x]), axis = 1)]

masked = []
for i in range(len(beds)):
    masked.append([])

for i in df["POS"]:
    for j in range(len(beds)):
        bed = beds[j]
        _ = bed[(bed["Start"] <= i) & (bed["End"] >= i)]
        if _.shape[0] == 1:
            _n = "_".join(_["Name"].values[0].split("_")[:-1])
            _ = bed[bed["Name"].str.contains(_n)].sort_values("Start")
            masked[j].extend(_.index.values) # Two primers per amplicon

for i in range(len(masked)):
    masked[i] = set(masked[i])

for i in range(len(out_paths)):
    txt = " ".join([str(j) for j in masked[i]])
    print(txt)
    with open(prefix+out_paths[i]+'_masked_primer_indices.txt', 'w') as f:
        f.write(txt+'\n')
