import pandas as pd
from pandarallel import pandarallel 
from pyliftover import LiftOver

pandarallel.initialize(progress_bar=True)

lo = LiftOver('hg19', 'hg38')

def conv(x, chrom_name, col_name):
    c = lo.convert_coordinate(x[chrom_name], x[col_name])
    if not c:
        return [pd.NA, pd.NA]
    return c[0][:2]

def conv_start1(x):
    return conv(x, "Chrom1", "Start1")

def conv_end1(x):
    return conv(x, "Chrom1", "End1")

def conv_start2(x):
    return conv(x, "Chrom2", "Start2")

def conv_end2(x):
    return conv(x, "Chrom2", "End2")


clusters = pd.read_csv("~/BioData/chromatin_loops/GSM1872886_GM12878_CTCF_PET_clusters.txt", names=["Chrom1", "Start1", "End1", "Chrom2", "Start2", "End2", "Score"], sep='\t')
print(clusters.head())
clusters[["Chrom1_hg38","Start1_hg38"]] = clusters.parallel_apply(conv_start1, axis=1, result_type="expand")
clusters[["Chrom1_hg38","End1_hg38"]] = clusters.parallel_apply(conv_end1, axis=1, result_type="expand")
clusters[["Chrom2_hg38","Start2_hg38"]] = clusters.parallel_apply(conv_start2, axis=1, result_type="expand")
clusters[["Chrom2_hg38","End2_hg38"]] = clusters.parallel_apply(conv_end2, axis=1, result_type="expand")

print(clusters.head())

print(len(clusters))
print(clusters["Start1_hg38"].isna().sum())
print(clusters["End1_hg38"].isna().sum())
print(clusters["Start2_hg38"].isna().sum())
print(clusters["End2_hg38"].isna().sum())
clusters = clusters.dropna()
print(len(clusters))

clusters = clusters[["Chrom1_hg38", "Start1_hg38", "End1_hg38", "Chrom2_hg38", "Start2_hg38", "End2_hg38"]]

clusters = clusters[clusters.Start1_hg38<=clusters.End1_hg38]
clusters = clusters[clusters.Start2_hg38<=clusters.End2_hg38]
print(len(clusters))

clusters.to_csv("~/BioData/chromatin_loops/GSM1872886_GM12878_CTCF_PET_clusters_hg38.txt", sep='\t', index=False, header=False)