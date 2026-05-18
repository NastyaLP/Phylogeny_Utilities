import sys, os
import time
import configparser
from src.common import read_input
from matplotlib.patches import Patch
import matplotlib.pyplot as plt

import logging


#############################################################
# adding fixed colour class or higher tax. level to the annotations
#############################################################

def parse_config(configpath):

    configfile = os.path.join(configpath)
    conf = configparser.ConfigParser()
    conf.optionxform = str
    conf.read(configfile)

    return conf

@read_input
def preselect_columns(anotfile, column, colorlist):

    counts = dict()

    head = next(anotfile)
    try:
        icolumn = head.index(column)
    except ValueError:
        try:
            icolumn = head.index(column).lower()
        except ValueError:
            logging.error(f"No custom column naimed '{column}' in the annotation, please "
                          f"change the column name for the argument --color_by_column")
            exit()
        pass
    for line in anotfile:

        counts[line[icolumn]] = counts.get(line[icolumn], 0) + 1

    i = 0

    fea2color = dict()

    for k, v in sorted(counts.items(), key=lambda item: item[1]):

        if i == len(colorlist):
            break
        fea2color[k] = colorlist[i]
        i += 1

    return fea2color


def plot_legend(fea2color, outdir, indir):

    fig, ax = plt.subplots()
    ax.axis('off')

    handles = [Patch(color=color, label=label) for label, color in fea2color.items()]

    # Add legend
    ax.legend(handles=handles, loc="center")
    if outdir:
        outpdf = os.path.join(outdir, "phylogeny_branch_color_legend.pdf")
    else:
        outdir = os.path.dirname(indir)
        logging.warning(f"No output directory provided. Legend plot 'phylogeny_branch_color_legend.pdf' can "
                     f" be found in {os.path.realpath(outdir)}")

        outpdf = os.path.join(outdir, "phylogeny_branch_color_legend.pdf")

    fig.savefig(outpdf)

def parse_anot(anotfile, colordict):

    tax_domain = "gtdb_domain"
    tax_class = "gtdb_domain"
    outhandle = open(os.path.splitext(anotfile)[0] + "_colormap_domain.tsv", "w")
    outcolor = open("../colours_tree.txt", "w")
    with open(anotfile) as inanot:

        header = next(inanot).strip("\n").split("\t")

        header.append("Colormap")
        outhandle.write("\t".join(header) + "\n")
        header = [el.lower() for el in header]
        idomain = header.index(tax_domain)
        iclass = header.index(tax_class)

        for line in inanot:

            line = line.strip("\n").split("\t")
            taxdomain = line[idomain]
            taxclass = line[iclass]

            if taxdomain == "Bacteria":
                color = colordict[taxdomain].replace('"', '')
                line.append(color)
            else:
                color = colordict[taxclass].replace('"', '')
                line.append(color)

            line = [el.replace(" ", "_") for el in line]

            outhandle.write("\t".join(line) + "\n")
            outcolor.write("\t".join([line[0], color]) + "\n")



