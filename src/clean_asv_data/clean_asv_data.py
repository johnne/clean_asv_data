#!/usr/bin/env python
import argparse
from argparse import ArgumentParser
import pandas as pd
import sys
import tqdm
from clean_asv_data.__main__ import read_config, generate_reader, read_clustfile


def read_counts(countsfile, blanks, chunksize=None, nrows=None):
    """
    Read the counts file in chunks, if list of blanks is given, count occurrence
    in blanks and return as a column <in_n_blanks>. Also calculate max and
    sum for each ASV.
    """
    reader = generate_reader(countsfile, chunksize=chunksize, nrows=nrows)
    sys.stderr.write("####\n" f"Reading counts from {countsfile}\n")
    dataframe = pd.DataFrame()
    for i, df in enumerate(tqdm.tqdm(reader,unit=" chunks",)):
        if i == 0:
            samples = df.shape[1]
        if len(blanks) > 0:
            blank_counts = df.loc[:, blanks]
            # calculate occurrence in blanks
            asv_blank_count = pd.DataFrame(
                blank_counts.gt(0).sum(axis=1), columns=["in_n_blanks"]
            )
            asv_blank_count["in_percent_blanks"] = (
                asv_blank_count.div(len(blanks)) * 100
            )
        # calculate ASV sum (remove blanks)
        asv_sum = pd.DataFrame(df.drop(blanks, axis=1).sum(axis=1), columns=["ASV_sum"])
        # calculate ASV max
        asv_max = pd.DataFrame(df.drop(blanks, axis=1).max(axis=1), columns=["ASV_max"])
        _dataframe = pd.merge(asv_sum, asv_max, left_index=True, right_index=True)
        _dataframe = pd.merge(
            _dataframe, asv_blank_count, left_index=True, right_index=True
        )
        dataframe = pd.concat([dataframe, _dataframe])
    sys.stderr.write(
        f"Read counts for {dataframe.shape[0]} ASVs in " f"{samples} samples\n"
    )
    return dataframe


def clean_by_taxonomy(dataframe, rank="Family"):
    """
    Removes ASVs if they are 'unassigned' at <rank> or <rank> contains '_X'
    """
    df = dataframe.copy()
    before = df.shape[0]
    sys.stderr.write("####\n" f"Removing ASVs unclassified at {rank}\n")
    cleaned = df.loc[
        (~df[rank].str.contains("_X+$")) & (~df[rank].str.startswith("unclassified"))
    ]
    after = cleaned.shape[0]
    sys.stderr.write(
        f"{before - after} ASVs removed, {cleaned.shape[0]} ASVs remaining\n"
    )
    return cleaned


def read_blanks(f):
    sys.stderr.write("####\n" f"Reading list of blanks from {f}\n")
    with open(f, "r") as fhin:
        blanks = [x.rstrip() for x in fhin.readlines()]
    sys.stderr.write(f"{len(blanks)} blanks read\n")
    return blanks


def clean_by_reads(dataframe, min_clust_count=3):
    """
    Remove clusters with a sum less than <min_reads> across samples
    """
    sys.stderr.write(
        "####\n" f"Removing ASVs in clusters with <{min_clust_count} total reads\n"
    )
    df = dataframe.copy()
    # Groupby cluster column and sum values in ASV_sum
    cl_sum = df.groupby("cluster").sum(numeric_only=True)
    # Get list of clusters to remove
    cl_to_remove = cl_sum.loc[cl_sum["ASV_sum"] < min_clust_count].index
    # Get list of ASVs in said clusters
    asvs_to_remove = df.loc[df["cluster"].isin(cl_to_remove)].index
    before = df.shape[0]
    df.drop(asvs_to_remove, inplace=True)
    after = df.shape[0]
    sys.stderr.write(f"{before - after} ASVs removed, {df.shape[0]} ASVs remaining\n")
    return df


def clean_by_blanks(dataframe, blanks, mode="asv", max_blank_occurrence=5):
    """
    Removes clusters with ASVs present in > <max_blank_occurrence>% of blanks
    """
    sys.stderr.write(
        "####\n" f"Removing {mode}s in >{max_blank_occurrence}% of blanks\n"
    )
    df = dataframe.copy()
    before = df.shape[0]
    # Get list of ASVs to remove
    to_remove = df.loc[df.in_percent_blanks > max_blank_occurrence].index
    if mode == "cluster":
        # Get clusters for said ASVs
        to_remove_cl = df.loc[to_remove, "cluster"]
        # Get all ASVs in said clusters
        to_remove = list(df.loc[df["cluster"].isin(list(to_remove_cl))].index)
    df.drop(to_remove, inplace=True)
    after = df.shape[0]
    sys.stderr.write(f"{before - after} ASVs removed, {df.shape[0]} ASVs remaining\n")
    return df


def main(args):
    # Read config
    args = read_config(args.configfile, args)
    print(args.nrows)
    # Read taxonomy + clusters
    asv_taxa = read_clustfile(args.clustfile)
    # Read blanks
    blanks = read_blanks(args.blanksfile) if args.blanksfile else []
    # Read counts
    counts = read_counts(
        countsfile=args.countsfile, blanks=blanks, chunksize=args.chunksize, nrows=args.nrows
    )
    # Clean by taxonomy
    asv_taxa_cleaned = clean_by_taxonomy(dataframe=asv_taxa, rank=args.clean_rank)
    # Merge counts + taxonomy
    asv_taxa_cleaned = pd.merge(
        asv_taxa_cleaned, counts, left_index=True, right_index=True
    )
    # Clean by blanks
    if args.blanksfile:
        asv_taxa_cleaned = clean_by_blanks(
            dataframe=asv_taxa_cleaned,
            blanks=args.blanksfile,
            mode=args.blank_removal_mode,
            max_blank_occurrence=args.max_blank_occurrence,
        )
    # Clean by read sum
    asv_taxa_cleaned = clean_by_reads(
        dataframe=asv_taxa_cleaned, min_clust_count=args.min_clust_count
    )
    # Write to output
    with sys.stdout as fhout:
        sys.stderr.write(
            "####\n" f"Writing {asv_taxa_cleaned.shape[0]} ASVs to stdout\n"
        )
        asv_taxa.index.name = "ASV"
        asv_taxa_cleaned.to_csv(fhout, sep="\t")


def main_cli():
    parser = ArgumentParser(
        """
        This script cleans clustering results by removing ASVs if:
        - unassigned or ambiguous taxonomic assignments (e.g.  
        'unclassified' or '_X' in rank labels) 
        - if belonging to clusters present in > max_blank_occurrence% of blanks
        - if belonging to clusters with < min_clust_count total reads
        """
    )
    io_group = parser.add_argument_group("input/output")
    io_group.add_argument("--countsfile", type=str, help="Counts file of ASVs")
    io_group.add_argument(
        "--clustfile",
        type=str,
        help="Taxonomy file for ASVs. Should also "
        "include a"
        "column with cluster designation.",
    )
    io_group.add_argument(
        "--blanksfile", type=str, help="File with samples that are 'blanks'"
    )
    io_group.add_argument("--output", type=str, help="Output file with cleaned results")
    params_group = parser.add_argument_group("params")
    params_group.add_argument("--configfile", type=str, default="config.yml",
                              help="Path to a yaml-format configuration file. Can be used to set arguments.")
    params_group.add_argument(
        "--clean_rank",
        type=str,
        help="Remove ASVs unassigned at this taxonomic " "rank (default Family)",
    )
    params_group.add_argument(
        "--max_blank_occurrence",
        type=int,
        help="Remove ASVs occurring in clusters where at "
        "least one member is present in "
        "<max_blank_occurrence>%% of blank samples. "
        "(default 5)",
    )
    params_group.add_argument(
        "--blank_removal_mode",
        type=str,
        choices=["cluster", "asv"],
        help="How to remove sequences based on "
        "occurrence in blanks. If 'asv' ("
        "default) remove "
        "only ASVs that occur in more than "
        "<max_blank_occurrence>%% of blanks. If "
        "'cluster', remove ASVs in clusters where "
        "one or more ASVs is above the "
        "<max_blank_occurrence> threshold",
    )
    params_group.add_argument(
        "--min_clust_count",
        type=int,
        help="Remove clusters with < <min_clust_count> "
        "summed across samples (default 3)",
    )
    debug_group = parser.add_argument_group("debug")
    debug_group.add_argument(
        "--chunksize",
        type=int,
        help="Size of chunks (in lines) to read from " "countsfile",
    )
    debug_group.add_argument(
        "--nrows",
        type=int,
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args()
    main(args)

