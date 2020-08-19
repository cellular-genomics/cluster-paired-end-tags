import argparse
import logging
import pandas as pd
import numpy as np
import numba
import functools
import itertools
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(message)s', datefmt='%H:%M:%S')


def add_build_args(parser):
    parser.add_argument("--pets_filename", type=str, action='append',
                        default=["~/BioData/4DNUCLEOME/4DNFI2BAXOSW_GM12878_CTCF_rep1_hiseq.bedpe"],
                        help=".bedpe files containing raw PETs (no header)")
    parser.add_argument("--clusters_filename", type=str,
                        default="~/BioData/chromatin_loops/4DNFI2BAXOSW_GM12878_CTCF_rep1_hiseq.bedpe.1.9.25.clusters",
                        help=".bedpe output file with clustered PETs (no header)")
    parser.add_argument("--self_ligation", type=int, default=8000, help="Self-ligation genomic span (default: 8000)")
    parser.add_argument("--extension", type=int, default=25,
                        help="No of base pairs to extend both ends of PETs (default: 25)")
    parser.add_argument("--pet_cutoff", type=int, default=1,
                        help="Minimum number of PET counts to take PET into consideration (default: 1)")
    parser.add_argument("--cluster_cutoff", type=int, default=9,
                        help="Minimum number of total counts to consider as a cluster (default: 9)")
    parser.add_argument("--nrows", type=int,
                        help="If provided limits the number of rows read from the PET file to nrows")


def cluster_PETs(args):
    # reading the file
    logging.info(f"Reading PETs from {args.pets_filename} ...")
    columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "cnt"]
    pets = pd.concat([pd.read_csv(f, sep="\t", header=None, names=columns, low_memory=False, nrows=args.nrows) for f in args.pets_filename])
    logging.info(f"Read {len(pets):,} PETs.")

    # pre-proccess
    logging.info(f"Preprocessing (Exstension: {args.extension}bp, Self-ligation genomic span: {args.self_ligation}bp, "
                 f"PET cutoff: {args.pet_cutoff}) ...")

    # check the data integrity
    invalid_pets = pets[(pets.chrom1 != pets.chrom2) | (pets.start1 > pets.end1) | (pets.start2 > pets.end2)]
    if len(invalid_pets) > 0:
        logging.info(f"{len(invalid_pets)} inter-chromosomal or misordered PETs are ignored:")
        logging.info(invalid_pets.head())
    pets = pets[(pets.chrom1 == pets.chrom2) & (pets.start1 <= pets.end1) & (pets.start2 <= pets.end2)]

    # keep only non-self-ligating PETs
    pets = pets[pets.start2 - pets.end1 >= args.self_ligation]

    # keep PETs with count >= PET cutoff
    pets = pets[pets.cnt >= args.pet_cutoff]

    # add extension
    pets.start1 = (pets.start1 - args.extension).clip(0)
    pets.end1 = pets.end1 + args.extension
    pets.start2 = (pets.start2 - args.extension).clip(0)
    pets.end2 = pets.end2 + args.extension

    logging.info("Done.")

    chroms = pets.groupby(["chrom1"]).size().to_dict()

    start1s = numba.typed.List([pets[pets.chrom1 == chrom].start1.to_numpy() for chrom in chroms.keys()])
    end1s = numba.typed.List([pets[pets.chrom1 == chrom].end1.to_numpy() for chrom in chroms.keys()])
    start2s = numba.typed.List([pets[pets.chrom1 == chrom].start2.to_numpy() for chrom in chroms.keys()])
    end2s = numba.typed.List([pets[pets.chrom1 == chrom].end2.to_numpy() for chrom in chroms.keys()])
    cnts = numba.typed.List([pets[pets.chrom1 == chrom].cnt.to_numpy() for chrom in chroms.keys()])

    step = 0
    changes = np.ones(shape=(len(chroms),))
    while np.sum(changes) > 0:
        # sorting# if changes[i] > 0 else None\
        logging.info(f"Sorting (step: #{step+1}, PETs: {len(pets):,}) ...")
        orders = numba.typed.List([np.lexsort((end2s[i], start2s[i], end1s[i], start1s[i]))
                                   for i in range(len(chroms))])
        logging.info("Done.")

        # clustering
        logging.info(f"Clustering (step: #{step+1}, PETs: {len(pets):,}) ...")

        @numba.jit(nopython=True, parallel=True)
        def cluster(chroms, old_changes, orders, start1s, end1s, start2s, end2s, cnts):
            changes = np.zeros(shape=(len(chroms)-1,), dtype=np.uint64)

            for idx in numba.prange(len(chroms)-1):
                if old_changes[idx] > 0:
                    order, start1, end1, start2, end2, cnt =\
                        orders[idx], start1s[idx], end1s[idx], start2s[idx], end2s[idx], cnts[idx]
                    for _i in range(chroms[idx]):
                        i = order[_i]
                        if cnt[i] == 0:
                            continue
                        _j = _i + 1
                        while _j < chroms[idx]:
                            j = order[_j]
                            if start1[j] > end1[i]:
                                break
                            if cnt[j] == 0:
                                _j += 1
                                continue
                            if ((start1[i] <= start1[j] and start1[j] <= end1[i]) or (start1[i] <= end1[j] and end1[j] <= end1[i])) and\
                               ((start2[i] <= start2[j] and start2[j] <= end2[i]) or (start2[i] <= end2[j] and end2[j] <= end2[i])):
                                start1[i] = min(start1[i], start1[j])
                                end1[i] = max(end1[i], end1[j])
                                start2[i] = min(start2[i], start2[j])
                                end2[i] = max(end2[i], end2[j])
                                cnt[i] += cnt[j]
                                cnt[j] = 0
                                changes[idx] += 1
                            _j += 1
            return changes
        changes = cluster(numba.typed.List(chroms.values()), changes, orders, start1s, end1s, start2s, end2s, cnts)
        logging.info(f"Done. Changes: {int(sum(changes)):,}")
        step += 1

    # save to file
    logging.info(f"Saving to {args.clusters_filename} (cluster cufoff: {args.cluster_cutoff})... ")
    pets = pd.DataFrame()
    for i, (chrom, size) in enumerate(chroms.items()):
        pets = pd.concat([
            pets, pd.DataFrame(data={"chrom1": itertools.repeat(chrom, size),
                                     "start1": start1s[i][orders[i]],
                                     "end1": end1s[i][orders[i]],
                                     "chrom2": itertools.repeat(chrom, size),
                                     "start2": start2s[i][orders[i]],
                                     "end2": end2s[i][orders[i]],
                                     "cnt": cnts[i][orders[i]]})])
    pets = pets[pets.cnt >= args.cluster_cutoff]
    pets.to_csv(args.clusters_filename, sep="\t", index=False, header=False)
    logging.info(f"Done. Saved {len(pets):,} clusters.")

    return pets


def main():
    parser = argparse.ArgumentParser()
    add_build_args(parser)
    args = parser.parse_args()
    cluster_PETs(args)


if __name__ == "__main__":
    main()
