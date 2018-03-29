import pandas as pd
import pysam

p = "/Users/karthik/hpc_downloads/2018.03.21/W113.trimmed.aligned.sorted.bam"
s = pysam.AlignmentFile(p, "rb")

cov = {
    "pos": [],
    "n": []
}

for pileupcolumn in s.pileup():
    cov["pos"].append(pileupcolumn.pos)
    cov["n"].append(pileupcolumn.n)
    for pileupread in pileupcolumn.pileups:
        if not pileupread.is_del and not pileupread.is_refskip:
            seq.append(pileupread.alignment.query_sequence[pileupread.query_position])

s.close()
with open("/Users/karthik/hpc_downloads/2018.03.21/W113.trimmed.aligned.sorted.fa", "w") as f:
    f.write(">W113.fa\n")
    f.write("".join(seq)+"\n")

cov = pd.DataFrame(cov)
cov["n"].plot(kind="bar")
plt.savefig("")
