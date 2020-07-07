import pandas as pd
import numpy as np

clusters = pd.read_csv("~/BioData/chromatin_loops/4DNFI2BAXOSW_GM12878_CTCF_rep1_hiseq.bedpe", sep="\t", names=["chrom1", "start1", "end1", "chrom2", "start2", "end2", "count"])
clusters = clusters[(clusters.chrom1=="chr8") & (clusters.start1>=57980000) & (clusters.end1<=59000000) & (clusters.start2>=57980000) & (clusters.end2<=59000000)]
clusters.start1 -= 100
clusters.start2 -= 100
clusters.end1 += 100
clusters.end2 += 100
clusters.to_csv("~/BioData/chromatin_loops/4DNFI2BAXOSW_GM12878_CTCF_rep1_hiseq_chr8_100.bedpe", sep="\t", header=False, index=False)
exit()

clusters = pd.read_csv("~/BioData/chromatin_loops/4DNFI2BAXOSW_GM12878_CTCF_hires.bedpe.2.4.clusters", sep="\t", names=["chrom1", "start1", "end1", "chrom2", "start2", "end2", "count"])

counts = clusters["count"].to_numpy()

print(counts.shape)

#import matplotlib.pyplot as plt
#plt.hist(counts, bins = 50)
#plt.show()

all = counts.sum()

h=100
hist = np.ndarray(shape=(h,))
for i in range(h):
    hist[i] = np.sum(np.where(counts==i, 1, 0))
    print(f"{i} {hist[i]} {100*hist[i]/all}%")