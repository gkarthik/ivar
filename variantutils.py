import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

def plot_variants_by_amplicon(vdfs, beds, filename, kind="scatter"):
    amplicon_uniq = beds[0]["Name"].apply(lambda x: "_".join(x.split("_")[:-1])).unique().tolist()
    d = int(np.round(len(amplicon_uniq)**0.5))
    f = plt.figure(figsize=(20,20))
    gs = gridspec.GridSpec(d, d)
    for _i, i in enumerate(amplicon_uniq):
        ax = plt.subplot(gs[_i])
        for _vdf, vdf in enumerate(vdfs):
            bed = beds[_vdf]
            _ = bed[bed["Name"].str.contains(i)].sort_values("Start")
            if _.shape[0] == 1:# To deal with trimmed off first and last primer.
                if _["Strand"].values[0] == "+":
                    vdf = vdf[vdf["POS"]>=_["Start"].values[0]]
                else:
                    vdf = vdf[vdf["POS"]<= _["End"].values[0]]
            else:
                coords = [_["Start"].values[0], _["End"].values[1]]
                vdf = vdf[(vdf["POS"]>=coords[0]) & (vdf["POS"]<=coords[1])]
            vdf["Percentage"] = (vdf["AD"]/vdf["DP"]) * 100
            if vdf.empty:
                continue
            if kind == "scatter":
                vdf.plot(x="POS", y="Percentage", ls="None", marker="o", ax = ax, alpha = 0.2, ms = 2)
            else:
                vdf["POS"].plot.hist(ax = ax, alpha = 0.5)
        ax.set_title(i)
        if ax.legend_ != None:
            ax.legend_.remove()
    plt.tight_layout()
    plt.savefig(filename)
    plt.clf()
    plt.close()
