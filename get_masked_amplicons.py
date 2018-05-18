#!/usr/local/bin/python3.5
import sys
import pandas as pd
from scipy.stats import fisher_exact as fe
from statsmodels.sandbox.stats.multicomp import multipletests

from variantutils import create_variant_dataframe, plot_variants_by_amplicon

prefix = sys.argv[1]
freq = float(sys.argv[2])
bed_path = sys.argv[3]
df_paths = sys.argv[4:]

bed = pd.read_table(bed_path, sep="\t", names=["Region", "Start", "End", "Name", "Score", "Strand"])

vdfs = []
for _i,i in enumerate(df_paths):
    _ = pd.read_table(i, comment = "#", compression="gzip", names = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE"])
    _ = create_variant_dataframe(_[(_["ALT"] != ".") & (_["ALT"]!=_["REF"])].reset_index(drop=True))
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
for i in df["POS"]:
    _ = bed[(bed["Start"] <= i) & (bed["End"] >= i)]
    if _.shape[0] == 1:
        _n = "_".join(_["Name"].values[0].split("_")[:-1])
        _ = bed[bed["Name"].str.contains(_n)].sort_values("Start")
        masked.extend(_.index.values) # Two primers per amplicon

masked = set(masked)
txt = " ".join([str(i) for i in masked])
print(txt)

with open(prefix+'_masked_primer_indices.txt', 'w') as f:
    f.write(txt+'\n')
