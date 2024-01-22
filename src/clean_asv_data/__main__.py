import pandas as pd
import yaml
import os
import importlib.resources
import sys


class objectview(object):
    def __init__(self, d):
        self.__dict__ = d


def testdata():
    """
    Generates dataframe with 10 ASVs in 3 clusters:

             cluster Phylum Class Order           Family            Genus Species
    ASV_ID
    ASV1    cluster1     PA    CA    OA               FA  unclassified.FA      SA
    ASV2    cluster1     PA    CA    OA               FA               GA      SA
    ASV3    cluster1     PA    CA    OA               FA               GA      SA
    ASV4    cluster2     PB    CB    OB  unclassified.OA               GB      SB
    ASV5    cluster2     PB    CB    OB             FB_X               GB      SB
    ASV6    cluster2     PB    CB    OB             FB_X               GB      SB
    ASV7    cluster3     PC    CC    OC            FC_XX               GC      SC
    ASV8    cluster3     PC    CC    OC               FC               GC      SC
    ASV9    cluster3     PC    CC    OC               FC               GC      SC
    ASV10   cluster3     PC    CC    OC               FC               GC      SC

    Cleaning the dataframe by taxonomy, removing ASVs unclassified and ambiguous
    at Family would result in the following dataframe:

             cluster Phylum Class Order           Family            Genus Species
    ASV_ID
    ASV1    cluster1     PA    CA    OA               FA  unclassified.FA      SA
    ASV2    cluster1     PA    CA    OA               FA               GA      SA
    ASV3    cluster1     PA    CA    OA               FA               GA      SA
    ASV8    cluster3     PC    CC    OC               FC               GC      SC
    ASV9    cluster3     PC    CC    OC               FC               GC      SC
    ASV10   cluster3     PC    CC    OC               FC               GC      SC

    Skipping cleaning of ambiguous ASVs at Family would result in the following
    dataframe:
             cluster Phylum Class Order           Family            Genus Species
    ASV_ID
    ASV1    cluster1     PA    CA    OA               FA  unclassified.FA      SA
    ASV2    cluster1     PA    CA    OA               FA               GA      SA
    ASV3    cluster1     PA    CA    OA               FA               GA      SA
    ASV5    cluster2     PB    CB    OB             FB_X               GB      SB
    ASV6    cluster2     PB    CB    OB             FB_X               GB      SB
    ASV7    cluster3     PC    CC    OC            FC_XX               GC      SC
    ASV8    cluster3     PC    CC    OC               FC               GC      SC
    ASV9    cluster3     PC    CC    OC               FC               GC      SC
    ASV10   cluster3     PC    CC    OC               FC               GC      SC

    Skipping cleaning of unclassified ASVs at Family would result in the
    following:
             cluster Phylum Class Order           Family            Genus Species
    ASV_ID
    ASV1    cluster1     PA    CA    OA               FA  unclassified.FA      SA
    ASV2    cluster1     PA    CA    OA               FA               GA      SA
    ASV3    cluster1     PA    CA    OA               FA               GA      SA
    ASV4    cluster2     PB    CB    OB  unclassified.OA               GB      SB
    ASV8    cluster3     PC    CC    OC               FC               GC      SC
    ASV9    cluster3     PC    CC    OC               FC               GC      SC
    ASV10   cluster3     PC    CC    OC               FC               GC      SC


    :return: dataframe with cluster membership and taxonomic assignments
    """
    df = pd.DataFrame(
        data={
            "cluster": [
                "cluster1",
                "cluster1",
                "cluster1",
                "cluster2",
                "cluster2",
                "cluster2",
                "cluster3",
                "cluster3",
                "cluster3",
                "cluster3",
            ],
            "Phylum": ["PA", "PA", "PA", "PB", "PB", "PB", "PC", "PC", "PC", "PC"],
            "Class": ["CA", "CA", "CA", "CB", "CB", "CB", "CC", "CC", "CC", "CC"],
            "Order": ["OA", "OA", "OA", "OB", "OB", "OB", "OC", "OC", "OC", "OC"],
            "Family": [
                "FA",
                "FA",
                "FA",
                "unclassified.OA",
                "FB_X",
                "FB_X",
                "FC_XX",
                "FC",
                "FC",
                "FC",
            ],
            "Genus": [
                "unclassified.FA",
                "GA",
                "GA",
                "GB",
                "GB",
                "GB",
                "GC",
                "GC",
                "GC",
                "GC",
            ],
            "Species": ["SA", "SA", "SA", "SB", "SB", "SB", "SC", "SC", "SC", "SC"],
        },
        index=[
            "ASV1",
            "ASV2",
            "ASV3",
            "ASV4",
            "ASV5",
            "ASV6",
            "ASV7",
            "ASV8",
            "ASV9",
            "ASV10",
        ],
    )
    df.index.name = "ASV_ID"
    return df.set


def update_args(args, config):
    """
    Overrides arguments set in configfile with arguments from commandline

    :param args:
    :param config:
    :return:
    """
    for key, value in args.__dict__.items():
        # If argument is not set in configfile, add to config dict
        if key not in config.keys():
            config[key] = value
        # If argument is set in configfile, update config dict if value is not None
        else:
            if value is not None:
                config[key] = value
    return config


def load_configfile(configfile):
    with open(configfile, "r") as fhin:
        return yaml.safe_load(fhin)


def read_config(configfile, args):
    """
    Stores parameters from a yaml configfile in a dictionary

    Any arg given on commandline overrides the config settings. The dictionary
    is turned back into an object with attributes just like the args.

    :param configfile: Path to a yaml config file
    :return: config dict
    """
    # Read default config from package
    default_configfile = str(
        importlib.resources.files("clean_asv_data") / "config/config.yml"
    )
    config = load_configfile(default_configfile)
    if os.path.exists(configfile):
        cl_config = load_configfile(configfile)
    else:
        cl_config = {}
    # Update default config with config from cmd
    config.update(cl_config)
    config = update_args(args, config)
    args = objectview(config)
    return args


def read_blanks(f=None):
    if f is None:
        return []
    sys.stderr.write("####\n" f"Reading list of blanks from {f}\n")
    with open(f, "r") as fhin:
        blanks = [x.rstrip() for x in fhin.readlines()]
    sys.stderr.write(f"{len(blanks)} blanks read\n")
    return blanks


def read_metadata(f, index_name="sampleID_SEQ"):
    sys.stderr.write("####\n" f"Reading metadata from {f}\n")
    df = pd.read_csv(f, sep="\t", header=0, comment="#")
    return df.set_index(index_name)


def read_clustfile(f, sep="\t"):
    """
    Reads a cluster membership file for ASVs

    :param f:
    :param sep:
    :return:
    """
    return pd.read_csv(f, sep=sep, index_col=0, header=0)


def generate_reader(f, chunksize, nrows):
    """
    Sets up a reader with pandas. Handles both chunksize>=1 and chunksize=None

    :param f: Input file
    :param chunksize: Number of rows to read per chunk
    :param nrows: Number of total rows to read
    :return:
    """
    if nrows == 0:
        nrows = None
    r = pd.read_csv(
        f, sep="\t", index_col=0, header=0, nrows=nrows, chunksize=chunksize
    )
    if chunksize is not None:
        return r
    return [r]
