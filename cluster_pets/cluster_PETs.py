import argparse
import logging
import pandas as pd
import numpy as np
import numba
import functools
logging.basicConfig(level=logging.INFO)


def add_build_args(parser):
    parser.add_argument("--pets_filename", type=str,
                        default="~/BioData/chromatin_loops/4DNFI2BAXOSW_GM12878_CTCF_rep1_hiseq.bedpe",
                        help=".bedpe file containing raw PETs (no header)")
    parser.add_argument("--clusters_filename", type=str,
                        default="~/BioData/chromatin_loops/4DNFI2BAXOSW_GM12878_CTCF_rep1_hiseq.bedpe.1.9.25.py3.clusters",
                        help=".bedpe output file with clustered PETs (no header)")
    parser.add_argument("--self_ligation", type=int, default=8000, help="Self-ligation genomic span (default: 8000)")
    parser.add_argument("--extension", type=int, default=25,
                        help="No of base pairs to extend both ends of PETs (default: 25)")
    parser.add_argument("--pet_cutoff", type=int, default=1,
                        help="Minimum number of PET counts to take PET into consideration (default: 1)")
    parser.add_argument("--cluster_cutoff", type=int, default=9,
                        help="Minimum number of total counts to consider as a cluster (default: 9)")


def cluster_PETs(args):
    # reading the file
    logging.info(f"Reading PETs from {args.pets_filename} ...")
    columns = ["chrom1", "start1", "end1", "chrom2", "start2", "end2", "cnt"]
    pets = pd.read_csv(args.pets_filename, sep="\t", header=None, names=columns, low_memory=False)
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

    chrom, start1, end1, start2, end2, cnt = \
        pets.chrom1.to_numpy(dtype=str, copy=True), pets.start1.to_numpy(copy=True), pets.end1.to_numpy(copy=True),\
        pets.start2.to_numpy(copy=True), pets.end2.to_numpy(copy=True), pets.cnt.to_numpy(copy=True)

    chrom_idxs = None
    step = 0
    while True:
        # sorting
        logging.info(f"Sorting (step: #{step+1}, PETs: {len(pets):,}) ...")

        @numba.jit(nopython=True)
        def calc_chrom_idxs(chrom, order):
            chrom_idxs = numba.typed.List()
            old_chrom = chrom[0]
            for _i in range(len(chrom)):
                i = order[_i]
                if _i == 0 or chrom[i] != old_chrom:
                    chrom_idxs.append(_i)
                    old_chrom = chrom[i]
            chrom_idxs.append(len(chrom))
            return chrom_idxs

        @numba.jit(nopython=True)
        def cmp(i, j):
            if chrom[i] != chrom[j]:
                return -1 if chrom[i] < chrom[j] else 1
            if start1[i] != start1[j]:
                return start1[i] - start1[j]
            if end1[i] != end1[j]:
                return end1[i] - end1[j]
            if start2[i] != start2[j]:
                return start2[i] - start2[j]
            if end2[i] != end2[j]:
                return end2[i] - end2[j]
            return 0

        order = list(range(len(chrom)))
        order.sort(key=functools.cmp_to_key(cmp))
        order = np.array(order)
        logging.info("Done.")

        if not chrom_idxs:
            chrom_idxs = calc_chrom_idxs(chrom, order)

        # clustering
        logging.info(f"Clustering (step: #{step+1}, PETs: {len(pets):,}) ...")

        @numba.jit(nopython=True, parallel=True)
        def cluster(chrom_idxs, order, start1, end1, start2, end2, cnt):
            changes = np.zeros(shape=(len(chrom_idxs)-1,), dtype=np.uint64)

            for idx in numba.prange(len(chrom_idxs)-1):
                for _i in range(chrom_idxs[idx], chrom_idxs[idx+1]):
                    i = order[_i]
                    if cnt[i] == 0:
                        continue
                    _j = _i + 1
                    while _j < chrom_idxs[idx+1]:
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
        changes = cluster(chrom_idxs, order, start1, end1, start2, end2, cnt)
        sum_changes = int(sum(changes))
        logging.info(f"Done. Changes: {sum_changes:,}")
        if sum_changes == 0:
            break
        step += 1

    # save to file
    logging.info(f"Saving to {args.clusters_filename} (cluster cufoff: {args.cluster_cutoff})... ")
    pets = pd.DataFrame(data={"chrom1": chrom[order], "start1": start1[order], "end1": end1[order],
                              "chrom2": chrom[order], "start2": start2[order], "end2": end2[order],
                              "cnt": cnt[order]})
    pets = pets[pets.cnt >= args.cluster_cutoff]
    pets.to_csv(args.clusters_filename, sep="\t", index=False, header=False)
    logging.info(f"Done. Saved {len(pets)} clusters.")

    return pets


def main():
    parser = argparse.ArgumentParser()
    add_build_args(parser)
    args = parser.parse_args()
    cluster_PETs(args)


if __name__ == "__main__":
    main()
