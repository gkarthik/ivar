import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

def create_variant_dataframe(df):
    # Make each variant a new row
    d = {
        "POS": [],
        "REF": [],
        "ALT": [],
        "DP": [],
        "AD": []
    }
    for _i, i in enumerate(df["ALT"].str.split(",")):
        ref = df.ix[_i]["REF"]
        d["ALT"].extend([j for j in i])
        d["REF"].extend([ref] * len(i))
        d["POS"].extend([df.ix[_i]["POS"]] * len(i))
        ad = [int(j) for j in df.ix[_i]["SAMPLE"].split(":")[-1].split(",")]
        d["AD"].extend(ad[1:])
        dp = sum(ad)
        d["DP"].extend([dp] * len(i))
    return pd.DataFrame(d)

def plot_variants_by_amplicon(vdfs, bed, filename, kind="scatter"):
    amplicon_uniq = bed["Name"].apply(lambda x: "_".join(x.split("_")[:-1])).unique().tolist()
    d = int(np.round(len(amplicon_uniq)**0.5))
    f = plt.figure(figsize=(20,20))
    gs = gridspec.GridSpec(d, d)
    for _i, i in enumerate(amplicon_uniq):
        ax = plt.subplot(gs[_i])
        _ = bed[bed["Name"].str.contains(i)].sort_values("Start")
        coords = [_["Start"].values[0], _["End"].values[1]]
        for vdf in vdfs:
            vdf = vdf[(vdf["POS"]>=coords[0]) & (vdf["POS"]<=coords[1])]
            vdf["Percentage"] = (vdf["AD"]/vdf["DP"]) * 100
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
