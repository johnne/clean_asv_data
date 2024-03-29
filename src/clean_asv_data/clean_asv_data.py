#!/usr/bin/env python
import argparse
from argparse import ArgumentParser
import pandas as pd
import sys
import os
import tqdm
from clean_asv_data.__main__ import (
    read_config,
    generate_reader,
    read_clustfile,
    read_metadata,
)


def read_counts(
    countsfile,
    metadata=None,
    split_col="dataset",
    blanks=None,
    chunksize=None,
    nrows=None,
):
    """
    Read the counts file in chunks, if list of blanks is given, count occurrence
    in blanks and return as a column <in_n_blanks>. Also calculate max and
    sum for each ASV.
    """
    if blanks is None:
        blanks = []
    reader = generate_reader(countsfile, chunksize=chunksize, nrows=nrows)
    sys.stderr.write("####\n" f"Reading counts from {countsfile}\n")
    data = {}
    warnings = []
    n_asvs = n_samples = 0
    n_datasets = []
    for i, df in enumerate(
        tqdm.tqdm(
            reader,
            unit=" chunks",
        )
    ):
        n_asvs += df.shape[0]
        if i == 0:
            n_samples = df.shape[1]
            sample_names = list(df.columns)
            if metadata is None:
                # set up dummy metadata
                split_col = "dataset"
                metadata = pd.DataFrame(
                    data={split_col: [split_col] * len(sample_names)},
                    index=sample_names,
                )
        # get unique values of split_col
        split_col_vals = metadata[split_col].unique()
        for val in split_col_vals:
            if val not in data.keys():
                data[val] = pd.DataFrame()
            # get samples corresponding to split
            val_samples = metadata.loc[metadata[split_col] == val].index
            # get intersection of val_samples and the df columns
            val_samples_intersect = list(set(val_samples).intersection(df.columns))
            if len(val_samples_intersect) == 0 and i == 0:
                warnings.append(
                    "####\n" f"No samples found in counts data for {val}, skipping...\n"
                )
                continue
            n_datasets.append(val)
            # get samples in val_samples missing from df columns
            missing_samples = set(val_samples).difference(val_samples_intersect)
            if len(missing_samples) > 0 and i == 0:
                warnings.append(
                    "####\n"
                    f"WARNING: {len(missing_samples)} samples in metadata file are missing from counts file for {val}:\n"
                )
                warnings.append(f"{', '.join(missing_samples)} \n")
            # split the counts dataframe
            split_df = df.loc[:, val_samples_intersect]
            # get blanks for this dataset
            val_blanks = list(set(blanks).intersection(val_samples_intersect))
            # get blanks in metadata actually present in counts data
            val_blanks_intersect = list(
                set(val_blanks).intersection(val_samples_intersect)
            )
            # get blanks missing
            # missing_blanks = set(val_blanks).difference(val_blanks_intersect)
            if i == 0:
                warnings.append(
                    "####\n" f"{len(val_blanks_intersect)} blanks found for {val}\n"
                )
            # calculate ASV sum (remove blanks)
            asv_sum = pd.DataFrame(
                split_df.drop(val_blanks_intersect, axis=1).sum(axis=1),
                columns=["ASV_sum"],
            )
            # calculate ASV max
            asv_max = pd.DataFrame(
                split_df.drop(val_blanks_intersect, axis=1).max(axis=1),
                columns=["ASV_max"],
            )
            _dataframe = pd.merge(asv_sum, asv_max, left_index=True, right_index=True)
            if len(val_blanks_intersect) > 0:
                blank_counts = split_df.loc[:, val_blanks_intersect]
                # calculate occurrence in blanks
                asv_blank_count = pd.DataFrame(
                    blank_counts.gt(0).sum(axis=1), columns=["in_n_blanks"]
                )
                asv_blank_count["in_percent_blanks"] = (
                    asv_blank_count.div(len(val_blanks_intersect)) * 100
                )
                _dataframe = pd.merge(
                    _dataframe, asv_blank_count, left_index=True, right_index=True
                )
            data[val] = pd.concat([data[val], _dataframe])
    sys.stderr.write(
        f"Read counts for {n_asvs} ASVs in "
        f"{n_samples} samples and {len(set(n_datasets))} datasets\n"
    )
    for item in warnings:
        sys.stderr.write(item)
    return data


def clean_by_taxonomy(dataframe, skip_ambig=False, skip_unclass=False, rank="Family"):
    """
    Removes ASVs if they are 'unassigned' at <rank> or <rank> contains '_X'
    """
    df = dataframe.copy()
    before = df.shape[0]
    cl_before = len(df["cluster"].unique())
    if skip_ambig and skip_unclass:
        sys.stderr.write("####\n" "Skipping cleaning by taxonomy\n")
        return df
    if not skip_ambig:
        sys.stderr.write("####\n" f"Removing ASVs ambiguous at {rank}\n")
        cleaned = df.loc[~df[rank].str.contains("_X+$")]
        n_ambig = before - cleaned.shape[0]
        cl_ambig = cl_before - len(cleaned["cluster"].unique())
        sys.stderr.write(f"{n_ambig} ASVs removed ({cl_ambig} clusters)\n")
    else:
        n_ambig = 0
        cl_ambig = 0
        cleaned = df.copy()
    if not skip_unclass:
        sys.stderr.write("####\n" f"Removing ASVs unclassified at {rank}\n")
        cleaned = cleaned.loc[~cleaned[rank].str.startswith("unclassified")]
        n_unclass = before - cleaned.shape[0] - n_ambig
        cl_unclass = cl_before - len(cleaned["cluster"].unique()) - cl_ambig
        sys.stderr.write(f"{n_unclass} ASVs removed ({cl_unclass} clusters)\n")
    else:
        n_unclass = 0
        cl_unclass = 0

    after = cleaned.shape[0]
    cl_after = len(cleaned["cluster"].unique())
    sys.stderr.write(
        f"{before - after} ASVs removed, {cleaned.shape[0]} ASVs remaining\n"
    )
    sys.stderr.write(
        f"{cl_before-cl_after} clusters removed, {cl_after} clusters remaining\n"
    )
    return cleaned


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


def clean_by_blanks(dataframe, blanks=None, mode="asv", max_blank_occurrence=5):
    """
    Removes ASVs present in > <max_blank_occurrence>% of blanks
    """
    if (
        blanks is None
        or len(blanks) == 0
        or "in_percent_blanks" not in dataframe.columns
    ):
        return dataframe
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
    data = {}
    # Read config
    args = read_config(args.configfile, args)
    if not args.output:
        outdir = "."
        output = "cleaned.tsv"
    else:
        outdir = os.path.dirname(args.output)
        if outdir=="":
            outdir = "."
        output = os.path.basename(args.output)
    # Read taxonomy + clusters
    sys.stderr.write("####\n" f"Reading clustfile {args.clustfile}\n")
    asv_taxa = read_clustfile(args.clustfile)
    sys.stderr.write(
        "###\n"
        f"Found {asv_taxa.shape[0]} ASVs in {len(asv_taxa['cluster'].unique())} clusters\n"
    )
    # Read metadata
    metadata = None
    blanks = None
    if args.metadata:
        metadata = read_metadata(args.metadata, index_name=args.metadata_index_name)
        # Extract blanks from metadata
        if not args.noblanks:
            blanks = list(
                metadata.loc[metadata[args.sample_type_col].isin(args.blank_val)].index
            )
            sys.stderr.write("####\n" f"Found {len(blanks)} blanks in metadata\n")
        else:
            blanks = []
    # Read counts (returns a dictionary)
    counts = read_counts(
        countsfile=args.countsfile,
        metadata=metadata,
        split_col=args.split_col,
        blanks=blanks,
        chunksize=args.chunksize,
        nrows=args.nrows,
    )
    # Clean by taxonomy
    asv_taxa_cleaned = clean_by_taxonomy(dataframe=asv_taxa, skip_ambig=args.skip_ambig, skip_unclass=args.skip_unclass, rank=args.clean_rank)
    # Merge counts + taxonomy
    for dataset, dataframe in counts.items():
        sys.stderr.write("####\n" f"Cleaning {dataset}\n")
        dataframe = pd.merge(
            asv_taxa_cleaned, dataframe, left_index=True, right_index=True
        )
        # Clean by blanks
        dataframe = clean_by_blanks(
            dataframe=dataframe,
            blanks=blanks,
            mode=args.blank_removal_mode,
            max_blank_occurrence=args.max_blank_occurrence,
        )
        # Clean by read sum
        dataframe = clean_by_reads(
            dataframe=dataframe, min_clust_count=args.min_clust_count
        )
        dataframe.index.name = "ASV"
        # Write to output
        if len(counts.keys()) > 1:
            outfile = f"{outdir}/{dataset}.{output}"
        else:
            outfile = f"{outdir}/{output}"
        with open(outfile, "w") as fhout:
            sys.stderr.write(
                "####\n"
                f"Writing cleaned {dataset} with {dataframe.shape[0]} ASVs to {outfile}\n"
            )
            dataframe.to_csv(fhout, sep="\t")


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
        "--metadata", type=str, help="Metadata file for splitting samples by datasets"
    )
    io_group.add_argument(
        "--metadata_index_name",
        type=str,
        help="Name of column in metadata file that contains sample ids",
        default="sampleID_NGI",
    )
    io_group.add_argument(
        "--split_col",
        type=str,
        help="Name of column in metadata file by which to split samples by prior to cleaning "
        "by blanks",
        default="dataset",
    )
    io_group.add_argument(
        "--output",
        type=str,
        help="Output file with cleaned results. If input data will be split "
        "into multiple datasets there will be one cleaned file per dataset"
        "with this parameter used to set the file name ending",
    )
    params_group = parser.add_argument_group("params")
    params_group.add_argument(
        "--configfile",
        type=str,
        default="config.yml",
        help="Path to a yaml-format configuration file. Can be used to set arguments.",
    )
    params_group.add_argument(
        "--sample_type_col",
        type=str,
        default="lab_sample_type",
        help="Use this column in metadata to identify sample type (default 'lab_sample_type')",
    )
    params_group.add_argument(
        "--blank_val",
        type=str,
        nargs="+",
        default=["buffer_blank", "extraction_neg", "pcr_neg"],
        help="Values in <sample_type_col> that identify blanks (default 'buffer_blank', 'extraction_neg', 'pcr_neg')",
    )
    params_group.add_argument(
        "--noblanks",
        action="store_true",
        help="Ignore blanks",
    )
    params_group.add_argument(
        "--clean_rank",
        type=str,
        help="Remove ASVs unassigned at this taxonomic " "rank (default Family)",
    )
    params_group.add_argument(
        "--skip_ambig",
        action="store_true",
        help="Skip cleaning of ASVs ambiguous at <clean_rank>. Ambiguous assignments end in '_X'")
    params_group.add_argument(
        "--skip_unclass",
        action="store_true",
        help="Skip cleaning of ASVs unclassified at <clean_rank>",
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
        default="asv",
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
