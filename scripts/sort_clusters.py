import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Sort clusters')
parser.add_argument('input')
parser.add_argument('output')
args = parser.parse_args()

columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "count"]

clusters = pd.read_csv(args.input, sep="\t", names=columns)
clusters=clusters[(clusters.chrom1=="chr1")&(clusters.chrom2=="chr1")]
clusters = clusters.sort_values(by=columns)
clusters.to_csv(args.output, sep="\t", index=False, header=False)