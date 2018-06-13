#!/usr/local/bin/python3.5

import sys
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import fisher_exact as fe
import numpy as np
from statsmodels.sandbox.stats.multicomp import multipletests
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.colors import Normalize
import matplotlib.lines as mlines
from matplotlib.collections import LineCollection
import argparse
from variantutils import plot_variants_by_amplicon

prefix = sys.argv[1]
freq = float(sys.argv[2])
df_bed_paths = sys.argv[3:]
df_paths = df_bed_paths[0::2]
bed_paths = df_bed_paths[1::2]

beds = []
for i in bed_paths:
    beds.append(pd.read_table(i, sep="\t", names=["Region", "Start", "End", "Name", "Score", "Strand"]))

vdfs = []
for _i,i in enumerate(df_paths):
    _ = pd.read_table(i, names=["POS", "REF", "ALT", "AD", "REV", "DP", "QUAL"], skiprows = 1)
    vdfs.append(_)

plot_variants_by_amplicon(vdfs, beds, prefix+"_variants_by_amplicon.png")
plot_variants_by_amplicon(vdfs, beds, prefix+"_variants_by_amplicon_hist.png", kind="hist")

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

cnorm = Normalize(vmin=0, vmax=len(cols))
sm = cm.ScalarMappable(norm=cnorm, cmap=cm.plasma)

plt.figure(figsize=(10, 8))
gs = gridspec.GridSpec(4, 1, height_ratios = [0.2,1,1,1])

# Percentage Plots
ax1 = plt.subplot(gs[1])
cols = df.columns[df.columns.str.contains("AD")]
for _i, i in enumerate(cols):
    _ = 100 * df[i]/df[i.replace("AD", "DP")]
    df[i.replace("AD", "Percentage")] = _
    c = sm.to_rgba(_i)
    s = df[i.replace("AD",col+"_pval")] <= pval_threshold
    c = [c if i else (222/255,222/255,222/255,1) for i in s]
    ax1.scatter(df["POS"], _, marker="o", color = c, edgecolor = "#000000", zorder = 2, alpha = 0.8)
    p1 = np.vstack((df["POS"], [-5] * len(df["POS"]))).T
    p2 = np.vstack((df["POS"], _)).T
    p = np.column_stack((p1, p2)).reshape(p1.shape[0],2,2)
    line_segments = LineCollection(p, lw=1, color="#585858", linestyle='dashed', zorder = 1)
    ax1.add_collection(line_segments)

_ = df[df.columns[df.columns.str.contains("Percentage")]].max(axis = 1).index
for j in range(0, len(_)):
    ax1.text(df["POS"][_[j]], _[j], str(df["REF"][_[j]])+" -> "+str(df["ALT"][_[j]]), ha="left", va="bottom", zorder = 3)

ax1.axhline(threshold, color="red")
ax1.set_ylim([-5, 105])

bed = beds[0]
ax1.set_xlim([bed["Start"].min(), bed["End"].max()])

ax1.set_title("Variant Frequency vs Position")
ax1.set_xlabel("Position")
ax1.set_ylabel("Frequency in %")

ax2 = plt.subplot(gs[2])
pval_cols = df.columns[df.columns.str.contains(r"\b"+col+"\b*_pval")]
for _i, i in enumerate(pval_cols):
    _ = -1 * df[i]
    c = sm.to_rgba(_i)
    s = df[i.replace("_pval", "_oddsratio")] * 10
    ax2.scatter(df["POS"], _, marker="o", color = c, edgecolor = "#000000", zorder = 2, alpha = 0.6, s= s)
    p1 = np.vstack((df["POS"], [-1.1] * len(df["POS"]))).T
    p2 = np.vstack((df["POS"], _)).T
    p = np.column_stack((p1, p2)).reshape(p1.shape[0],2,2)
    line_segments = LineCollection(p, lw=1, color="#585858", linestyle='dashed', zorder = 1)
    ax2.add_collection(line_segments)

_ = df[pval_cols].min(axis = 1).index
for j in range(0, len(_)):
    ax2.text(df["POS"][_[j]], _[j], str(-1 * df["REF"][_[j]])+" -> "+str(df["ALT"][_[j]]), ha="left", va="bottom", zorder = 3)

ax2.axhline(-1 * pval_threshold, color="red")
ax2.set_ylim([-1.1, 0.1])
ax2.set_xlim([bed["Start"].min(), bed["End"].max()])
ax2.set_title("P-value of Fisher's test vs Position")
ax2.set_xlabel("Position")
ax2.set_ylabel("-1 * pvalue")

ax3 = plt.subplot(gs[3])
binwidth = 0.01
for i in pval_cols:
    ax3.hist(df[i], bins=np.arange(0, 1 + binwidth, binwidth), alpha = 0.5)

ax3.axvline(0.05, color="red")
ax3.set_title("Histogram of pvalues")
ax3.set_xlabel("P value")
ax3.set_ylabel("Count")

ax4 = plt.subplot(gs[0])
y = True
for i in bed[bed["Strand"] == "+"].sort_values("Start").index:
    n = "_".join(bed.ix[i]["Name"].split("_")[:-1])
    y = not y
    _ = bed[bed["Name"].str.contains(n)]
    if _.shape[0] == 1:         # To deal with trimmed off first primer
        continue
    else if _.shape[0] == 2:
        f = [_[_["Strand"] == "+"]["Start"].values[0], _[_["Strand"] == "+"]["End"].values[0]]
        r = [_[_["Strand"] == "-"]["Start"].values[0], _[_["Strand"] == "-"]["End"].values[0]]
        c = "black"
        ax4.plot([f[0], r[1]], [int(y),int(y)], color=c, lw = 2)
        ax4.plot([f[0], f[1]], [int(y),int(y)], color="steelblue", lw = 4)
        ax4.plot([r[0], r[1]], [int(y),int(y)], color="indianred", lw = 4)

ax4.yaxis.set_visible(False)
ax4.spines['top'].set_visible(False)
ax4.spines['right'].set_visible(False)
ax4.spines['bottom'].set_visible(False)
ax4.spines['left'].set_visible(False)

ax4.set_title("Amplicon Scheme")
ax4.set_xlim([bed["Start"].min(), bed["End"].max()])

df.to_csv(prefix+".csv")

# plt.suptitle("Threshold of "+str(threshold))
plt.tight_layout()
plt.savefig(prefix+"_report.pdf")
plt.clf()
plt.close()
