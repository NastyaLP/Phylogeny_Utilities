from Bio import SeqIO
import csv
from functools import wraps
import gzip
import bz2
import sys, os
import yaml

csv.field_size_limit(sys.maxsize)

def read_input(func):

    @wraps(func)
    def func_wrapper(thefile, *args, **kwargs):

        if thefile.endswith("gz"):
            thefile = gzip.open(thefile, 'rt')
        elif thefile.endswith("bz2"):
            thefile = bz2.open(thefile, 'rt')
        else:
            thefile = open(thefile, 'r')

        if "nocsv" in args:
            return func(thefile, *args, **kwargs)
        elif "fasta" in args:
            return func(SeqIO.parse(thefile, "fasta"), *args, **kwargs)
        else:
            return func(csv.reader(thefile, delimiter=kwargs.get("delimiter", "\t")), *args, **kwargs)

    return func_wrapper


def parse_config_yaml(configname):
    """
    extracts information which proteins to be considered present by phylogeny
    and which ones by cluster

    """
    dirname = os.path.dirname(os.path.abspath(__file__))

    configpath = os.path.join(dirname, configname)

    with open(configpath) as stream:
        try:
            configload = yaml.safe_load(stream)

        except yaml.YAMLError as exc:
            print(exc)

    return configload