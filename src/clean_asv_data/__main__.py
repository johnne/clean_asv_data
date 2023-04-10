import pandas as pd
import yaml
import os
import importlib.resources
import sys


class objectview(object):
    def __init__(self, d):
        self.__dict__ = d


def update_args(args, config):
    """
    Overrides arguments set in configfile with arguments from commandline

    :param args:
    :param config:
    :return:
    """
    for key, value in args.__dict__.items():
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
