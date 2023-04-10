import argparse
from argparse import ArgumentParser
import pandas as pd
from clean_asv_data.__main__ import generate_reader, read_clustfile, read_config
from tqdm import tqdm
import sys


def find_consensus_taxonomies(
    clustdf, clust_column, ranks, consensus_ranks, consensus_threshold
):
    cluster_taxonomies = {}
    cons_ranks_reversed = consensus_ranks.copy()
    cons_ranks_reversed.reverse()
    for cluster in tqdm(
        clustdf[clust_column].unique(),
        desc="finding consensus taxonomies",
        unit=" clusters",
    ):
        rows = clustdf.loc[clustdf[clust_column] == cluster]
        for rank in cons_ranks_reversed:
            # Sum ASV sums up to rank
            rank_sums = rows.groupby(rank).sum(numeric_only=True)
            # Calculate percent for rank labels
            rank_sums_percent = rank_sums.div(rank_sums.sum()) * 100
            # Sort percent values in descending order
            rank_sums_percent.sort_values("ASV_sum", ascending=False, inplace=True)
            # Create a list <above_thresh> of number of labels above threshold
            above_thresh = list(
                rank_sums_percent.loc[
                    rank_sums_percent["ASV_sum"] > consensus_threshold
                ].index
            )
            # If only one assignment is above threshold, use this lineage to resolve taxonomy
            if len(above_thresh) == 1:
                taxlabel = above_thresh[0]
                lineage = (
                    rows.loc[rows[rank] == above_thresh[0]][ranks]
                    .head(1)
                    .to_dict(orient="index")
                )
                break
        cluster_taxonomies[cluster] = list(lineage.values())[0]
        ranks_below = cons_ranks_reversed[0 : cons_ranks_reversed.index(rank)]
        for r in ranks_below:
            cluster_taxonomies[cluster][r] = f"unresolved.{taxlabel}"
    return pd.DataFrame(cluster_taxonomies).T


def sum_asvs(countsfile, chunksize=None, nrows=None):
    reader = generate_reader(f=countsfile, chunksize=chunksize, nrows=nrows)
    asv_sum = pd.DataFrame()
    for df in tqdm.tqdm(reader, unit="chunks"):
        _asv_sum = pd.DataFrame(df.sum(numeric_only=True, axis=1), columns=["ASV_sum"])
        asv_sum = pd.concat([asv_sum, _asv_sum])
    return asv_sum.sort_values(ascending=False)


def main(args):
    args = read_config(args.configfile, args)
    asv_sum = sum_asvs(
        countsfile=args.countsfile, chunksize=args.chunksize, nrows=args.nrows
    )
    clustdf = read_clustfile(args.clustfile, sep="\t")
    clustdf = clustdf.loc[:, [args.clust_column] + args.ranks]
    clustdf = pd.merge(asv_sum, clustdf, left_index=True, right_index=True)
    resolved = find_consensus_taxonomies(
        clustdf=clustdf,
        clust_column=args.clust_column,
        ranks=args.ranks,
        consensus_ranks=args.consensus_ranks,
        consensus_threshold=args.consensus_threshold,
    )
    with sys.stdout as fhout:
        resolved.to_csv(fhout, sep="\t")


def main_cli():
    parser = ArgumentParser()
    parser.add_argument("--countsfile", type=str, help="Counts file of ASVs")
    parser.add_argument(
        "--clustfile",
        type=str,
        help="Taxonomy file for ASVs. Should also include a column with cluster designation.",
    )
    parser.add_argument(
        "--configfile",
        type=str,
        default="config.yml",
        help="Path to a yaml-format configuration file. Can be used to set arguments.",
    )
    parser.add_argument("--ranks", nargs="+", help="Ranks to include in the output.")
    parser.add_argument(
        "--clust_column",
        type=str,
        help="Name of cluster column, e.g. 'cluster'",
    )
    parser.add_argument(
        "--consensus_threshold",
        type=int,
        help="Threshold (in %%) at which to assign taxonomy to a cluster",
    )
    parser.add_argument(
        "--consensus_ranks",
        nargs="+",
        help="Ranks to use for calculating consensus. Must be present in the clustfile.",
    )
    parser.add_argument(
        "--chunksize",
        type=int,
        help="If countsfile is very large, specify chunksize to read it in a number of lines at a time",
    )
    parser.add_argument("--nrows", type=int, help=argparse.SUPPRESS)
    args = parser.parse_args()
    main(args)
