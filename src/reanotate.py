import shutil
from src.nexus_to_nexusfigtree import SeqidAnot, parse_figtree_config


def parse_nexus_taxlabels(treefile, seqidanot):

    with open(treefile) as tree:

        format = next(tree)

        if format.strip() != "#NEXUS":
            raise Exception(f"The given file {treefile} is not in nexus format. "
                            f"Check your input directory and extension given for reannotation")
        outpath = "tmp.nexus"
        outfile = open(outpath, "w")
        outfile.write(format)
        start = False
        anotbool = True
        for line in tree:

            if line.startswith(";"):
                start = False

            if start:
                seqid = line.split("[")[0]
                seqidbase = seqid.strip().replace("'", "")
                curseqid, lineadd = seqidanot.func(seqidbase, anotbool)
                outfile.write(curseqid)

            if line.strip().startswith("taxlabels"):
                start = True
                outfile.write(line)
            """
            if line.strip().startswith("begin figtree"):
                break
            """
            if not start:
                outfile.write(line)
        """
        for line in attrlines:
            outfile.write(line)
        """

        outfile.flush()
        outfile.close()
        shutil.copy(outpath, treefile)








