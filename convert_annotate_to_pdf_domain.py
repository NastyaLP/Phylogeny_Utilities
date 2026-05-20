#!/usr/bin/env python
import logging
logging.basicConfig(level=logging.WARNING)
################################################################################################
""" convert_annotate_to_pdf.py: starts with iqtree output, converts to basic nexus
    annotates and colors the tips, annotates the phylogeny,converts to pdf """
###############################################################################################

import sys, os
import argparse
import subprocess

from src.pdf_util import pdfmerge, pdf_anotate_single
from src.nexus_to_nexusfigtree import SeqidAnot, parse_figtree_config, parse_nexus_newick, write_nexus_figtree
from src.collapse import reanotate_nexus
from src.steady_colors import parse_anot, parse_config, preselect_columns, plot_legend
from src.reanotate import parse_nexus_taxlabels


__author__ = "Anastasiia Padalko"
__license__ = "GPL"
__version__ = "1.0"


pdfpathlist = []
tree2anot2count = dict()


def find_height_width(ntax):

    if ntax > 2000:
        return '2200', '1800'

    elif ntax > 1000:
        return '1800', '1400'

    else:
        return '1200', '1000'


def figtree_pdf_convert(treepath, ntips):

    width, height = find_height_width(ntips)

    pdfpath = os.path.splitext(treepath)[0] + ".pdf"
    figpath = f"{os.path.dirname(os.path.abspath(__file__))}/figtree.jar"

    subprocess.check_output(["java", "-jar", figpath, "-graphic", "PDF", "-width",  width, "-height", height, treepath, pdfpath])

    return pdfpath

def figtree_convert_pdf_static_width(treepath):

    width = '1000'
    height = '1000'

    pdfpath = os.path.splitext(treepath)[0] + ".pdf"
    figpath = f"{os.path.dirname(os.path.abspath(__file__))}/figtree.jar"

    subprocess.check_output(
        ["java", "-jar", figpath, "-graphic", "PDF", "-width", width, "-height", height, treepath, pdfpath])

    return pdfpath

def parse_tree_dir(indir,  queryext, seqidanot, convert, outdir, splitbool, merge_annotations, reanotate, collapse,
                   taxonomy_to_collapse):

    if not outdir:
        logging.warning(f"Output directory is not provided. Processed phylogenies can be found in the "
                     f"can be found in the same directory with the input ones with the disered extension '{queryext}'")

    taxlines, attrlines = parse_figtree_config("nexusfigtree_set.txt")


    for root, dirs, files in os.walk(indir):
        for file in files:

            if file.endswith(queryext):

                seqidanot.seqid2subindex = dict()

                if collapse:
                    outfilepath = reanotate_nexus(seqidanot, os.path.join(root, file), reanotate, taxonomy_to_collapse)
                    if convert:
                        pdfpath = figtree_convert_pdf_static_width(outfilepath)
                        pdfpathlist.append(pdfpath)
                elif not collapse and reanotate:

                    logging.warning("Re-annotating figtree-nexus files, the new nexus files will replace the old ones")
                    outfilepath = parse_nexus_taxlabels(os.path.join(root, file), seqidanot)

                else:
                    parse_tree(os.path.join(root, file), seqidanot, outdir, taxlines,
                                                attrlines, splitbool, convert,  merge_annotations)


def get_model(filepath):

    dirname = os.path.dirname(filepath)

    model = ""

    for file in os.listdir(dirname):

        if file.endswith("iqtree"):

            with open(os.path.join(dirname, file)) as inlog:

                for line in inlog:

                    if not model and line.lower().startswith("best-fit model"):
                        model = line.strip().split(": ")[1].split(" ")[0]

    return model


def process_annotation_note(anot2count,  ntax, treemodel, treepath):

    freq = 0.05

    anotlines = []

    anotlines.append(f"Filename: {treepath}")
    if not treemodel:
        anotlines.append(f"treemodel: None\n")
    else:
        anotlines.append(f"treemodel: {treemodel}\n")

    for key, typedict in anot2count.items():

        typelist = []

        if len(typedict) > 10:

            for k, value in sorted(typedict.items(), key=lambda item: item[1], reverse=True):

                if value/ntax >= freq:

                    typelist.append(":".join([k, str(value)]))

            typelist = ",".join(typelist)



        else:
            typelist = ",".join([":".join([k, str(value)]) for k, value in typedict.items()])

        anotlines.append("\t".join([key, typelist]) + "\n")
    return anotlines


def write_nexus_figtree_pdf(outre,  lines, taxlines, attrlines, outfile, convert, seqidanot, treemodel,
                            merge_annotations):

    outhandle = open(outfile, 'w')
    ntax = len(outre)

    write_nexus_figtree(outre, lines, taxlines, attrlines, outhandle)
    outhandle.close()
    anotlines = process_annotation_note(seqidanot.anot2count, ntax, treemodel, outfile)
    if seqidanot.merge:

        treebase = os.path.basename(os.path.splitext(outfile)[0])
        tree2anot2count[treebase] = anotlines

    if convert:
        pdfpath = figtree_pdf_convert(outfile,  ntax)
        pdfpathlist.append(pdfpath)
        if not merge_annotations:
            pdf_anotate_single(pdfpath, anotlines)



def parse_tree(file, seqidanot, outdir, taxlines, attrlines, splitbool, convert, merge_annotations):

    treemodel = get_model(file)

    outre, lines = parse_nexus_newick(file, seqidanot, splitbool)

    outenzyme = ""

    if seqidanot.enzbool:
        outenzyme = "_" + sorted(seqidanot.enzyme2count.items(), key=lambda item: item[1], reverse=True)[0][0]
        outenzyme = outenzyme.replace("/", "_")
    if outdir:
        basefile = os.path.splitext(os.path.basename(file))[0]

        dirname = file.split("/")[-2]
        basefilelist = basefile.split(".")

        if basefilelist[0] in dirname:
            basefilelist[0] = dirname
            basefile = ".".join(basefilelist)

        outfile = os.path.join(outdir, basefile + outenzyme)

    else:
        outfile = os.path.splitext(file)[0] + outenzyme

    if splitbool:

        for ind, treelines in enumerate(lines):
            if ind > 0:
                curoutfile = outfile + f"_0{ind}_auto_domain.figtree"
            else:
                curoutfile = outfile + f"_auto_domain.figtree"

            write_nexus_figtree_pdf(outre, treelines,  taxlines, attrlines, curoutfile, convert, seqidanot, treemodel,
                                    merge_annotations)

    else:

        outfile = outfile + f"_auto_domain.figtree"

        write_nexus_figtree_pdf(outre, lines, taxlines, attrlines, outfile, convert, seqidanot, treemodel,
                                merge_annotations)

def strip_hex_quotes(colordict):

    newcolordict = dict()

    for k, v in colordict.items():

        newcolordict[k] = v.replace('"', '')

    return newcolordict


def main():

    parser = argparse.ArgumentParser(
        prog='convert_annotate_to_pdf_domain.py',
        description='starts with iqtree output, converts to basic nexus,anotates and colors the tips, converts to pdf,'
                    'basic usage:\n'
                    'python convert_anotate_to_pdf.py -a <annotation_file> --merge_pdf <input_dir> <treefile_extension>')

    parser.add_argument('indir', help="input directory of any depth with phylogeny files", type=str)
    parser.add_argument('extension', help="file extension of the phylogeny files you want to process")
    parser.add_argument('-o', '--outdir', help="output directory, default saved in the input directory", type=str)

    parser.add_argument('-a', '--annotation', help="annotation file, should be supplied with at least seqid "
                                                              "and colormap column with hexcodes")
    parser.add_argument('-c', '--convert', action='store_true', default=False, help="convert to pdf, default False")
    parser.add_argument('-f', '--figtree_attribures', help="file with the desired figtree attributes to use")

    parser.add_argument('-m', '--merge_pdf', action='store_true', default=False,
                        help='merge phylogeny pdf in a single text annotated document')

    parser.add_argument('--pdf_annotations',  help='annotation columns that will be used for summary note'
                                                    'on the pdf files,  default all attributes'
                                                    'with frequency above 0.05 per phylogeny members')

    parser.add_argument('--color_by_column', help='column name for phylogeny branches to be colored by attributes within column'
                                               'by default it uses "Gtdb_domain" if available to'
                                               'color phylogenies by Archaea or Bacteria affiliation')
    parser.add_argument('--color_config', help="if '.ini' config provided colors are selected according to the config mappings (see test/color_config_example) "
                                               "for the attributes in the desired/default column")


    parser.add_argument('--reanotate', action='store_true', default=False,
                        help="reanotating existing nexus/figtree files")
    parser.add_argument( "--taxonomy_to_collapse",
        choices=["Gtdb", "Ncbi"],  # only allow these two
        help="Choose taxonomy to collapse by. By default uses both. \n"
             "Columns/Annotations in figtree for Ncbi: \n"
             "'Ncbi_phylum', 'Ncbi_class', 'Ncbi_order', 'Ncbi_species'. "
             "Columns/Annotations in figtree for Gtdb: "
             "'Gtdb_phylum', 'Gtdb_class', 'Gtdb_order', 'Gtdb_species'")
    parser.add_argument("--collapse", action="store_true", default=False,
                        help="collapses cartooned clades derives taxonomic labels from order "
                             "to phylum based on clade homogeneity")
    parser.add_argument('--merge_annotations', action='store_true',  default=False,
                        help='put pdf annotations only on the merged pdf requires --merge-pdf to be used.'
                                'Default puts annotations on all the single pdfs')

    parser.add_argument('--split_trees', action='store_true', default=False, help='if there are multiple phylogenies'
                                             'per file, adding this flag will split them')
    args = parser.parse_args()

    defaultconfig = "src/colours_tree_gtdb_domain.ini"

    if args.color_by_column and args.color_by_column.lower() != "gtdb_domain":
        column = args.color_by_column

        if not args.color_config:
            conf = parse_config(defaultconfig)
            colorlist = [c.strip().replace('"', '')for c in conf["colorlist"]["palette"].splitlines() if c.strip()]

            colordict = preselect_columns(args.annotation, column, colorlist)
        else:
            colordict = strip_hex_quotes(dict(parse_config(args.color_config)["colormap"]))

    else:
        column = "Gtdb_domain"

        if args.color_config:

            colordict = strip_hex_quotes(dict(parse_config(args.color_config)["colormap"]))

        else:
            colordict = strip_hex_quotes(dict(parse_config(defaultconfig)["colormap"]))


    plot_legend(colordict, args.outdir, args.indir)

    if args.outdir:
        os.makedirs(args.outdir, exist_ok=True)

    if not args.annotation:
        logging.warning("No annotation file provided, following the figtree configuration")

    if args.merge_annotations and not args.merge_pdf:
        raise Exception("Use --merge_pdf flag with --merge_annotations for the specified"
                        "annotation columns in the merged pdf file")

    seqidanot = SeqidAnot(args.annotation, colordict, column)
    seqidanot.set_merge(args.merge_pdf, args.pdf_annotations)

    parse_tree_dir(args.indir, args.extension, seqidanot, args.convert, args.outdir,
                   args.split_trees, args.merge_annotations, args.reanotate, args.collapse, args.taxonomy_to_collapse)

    if args.merge_pdf:

        if args.outdir:
            outpath = os.path.join(args.outdir, "joined_phylogeny_domain.pdf")
        else:

            outdir = os.path.dirname(args.indir)
            logging.warning(f"Output directory is not provided. All phylogenies pdf file 'joined_phylogeny_domain.pdf' "
                            f"can be found in {os.path.realpath(outdir)}")

            outpath = os.path.join(outdir, "joined_phylogeny_domain.pdf")

        if args.merge_annotations:
            pdfmerge(tree2anot2count, outpath, True, pdfpathlist)
        else:
            pdfmerge(tree2anot2count, outpath, False, pdfpathlist)


if __name__ == '__main__':
    main()
