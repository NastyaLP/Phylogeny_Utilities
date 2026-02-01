import sys, os
from nexus_to_nexusfigtree import SeqidAnot, parse_figtree_config, get_round

curtaxonomy = None
def treeline_collapse(treeline, seqidanot, itree, basename, taxonomy):

    labelfile = open(f"{basename}_{itree}_labels.txt", "w")

    seqids = []
    labelfilelist = []
    newtreeline = ''
    seqidbool = False
    closedbool = False
    subindexbool = False
    curseqid = ""
    lastseqid = ""

    seqidanot.reset_merge()
    countroot = 0
    seqidcount = 0
    checkboot = False

    cartoonbool = ""
    collapseseqids = []
    listofcollapsed = []
    allcolapsed = []
    curly_open = "{"
    curly_close = "}"
    prevchar = ""

    global curtaxonomy
    if taxonomy == "Ncbi":
        speciesanot = "Ncbi_species"
        curtaxonomy = "Ncbi"
    elif taxonomy == "Gtdb":
        speciesanot = "Gtdb_species"
        curtaxonomy = "Gtdb"
    else:
        speciesanot = "Gtdb_species"

    for char in treeline:

        if char == "'":

            if not seqidbool:
                seqidbool = True
                newtreeline += char
                continue
            else:
                seqidbool = False
                seqidcount += 1

                seqids.append(curseqid)
                lastseqid = curseqid
                collapseseqids.append(curseqid)
                allcolapsed.append(curseqid)
                curseqid = ""
        elif char.isalpha() and not seqidbool and closedbool and not newtreeline[-1].isdigit():

                seqidbool = True
                newtreeline += f"'{char}"
                curseqid += char
                continue

        elif char == "[":
            closedbool = False
            if seqidbool:
                seqidbool = False
                seqids.append(curseqid)
                lastseqid = curseqid

                collapseseqids.append(curseqid)
                allcolapsed.append(curseqid)
                curseqid = ""
                newtreeline += f"'{char}"
            else:
                newtreeline += f"{char}"
            continue

        elif char == "&" and lastseqid:
            ntaxa2count, gtaxa2count, color, nodelabel, subindexbool = (
                count_labels([lastseqid], seqidanot, subindexbool, speciesanot))
            if nodelabel:
                newtreeline += f"&nodelabel={nodelabel},"
            else:
                newtreeline += f"&nodelabel=None,"


            continue

        elif char == "=":
            if newtreeline[-7:] == "cartoon":
                cartoonbool += char
            elif newtreeline[-7:] == "rounded":
                checkboot = True
            newtreeline += char
            continue
        elif char == "}":
            if cartoonbool:
                value = cartoonbool.split(",")[1].strip("}")
                countids = int(cartoonbool.split(",")[0].strip("={"))
                cartoonbool = ""
                i = 0
                printbool = False

                if len(collapseseqids) > countids:

                    collapseseqids = collapseseqids[(len(collapseseqids) - countids):]

                if not collapseseqids:
                    printbool = True
                    for el in reversed(listofcollapsed):
                        collapseseqids += el
                        i += 1

                        if len(collapseseqids) == countids:
                            break

                    listofcollapsed = listofcollapsed[:-i]
                    labelfilelist = labelfilelist[:-i]

                ntaxa2count, gtaxa2count, color, nodelabel, subindexbool = count_labels(collapseseqids, seqidanot,
                                                                                        subindexbool, speciesanot)
                listofcollapsed.append(collapseseqids)

                newtreeline += f'{char}, !collapse={curly_open}"collapsed", {value}{curly_close},!color={color}, node={nodelabel}'

                newtreeline = go_back_newtreeline(newtreeline, nodelabel, collapseseqids)

                if not newtreeline.strip().startswith("tree tree0 = [&R]"):

                    raise Exception("NEXUS tree structure is corrupted")

                labelfilelist.append(f"ncbi_{ntaxa2count}\tgtdb_{gtaxa2count}\n")

                collapseseqids = []

            continue

        prevchar = char
        if cartoonbool:
            cartoonbool += char

        elif char == "]" or char == ",":
            if char == "]":
                closedbool = True
            if checkboot:

                checkboot = False

                newtreeline = newtreeline[:-2] + str(get_round(newtreeline[-2:].strip("="))[0])

            newtreeline += char
            #newtreeline += "'"
            openbool = True
            continue

        if seqidbool:
            curseqid += char

        newtreeline += char
        if char == ";":
            newtreeline += "\n"
            break
    for line in labelfilelist:
        labelfile.write(line)
    labelfile.flush()
    labelfile.close()
    #print(newtreeline)
    return seqids,  newtreeline

def go_back_newtreeline(treeline, nodelabel, collseqid):

    seqidbool = False
    seqidcount = 0
    newtreeline = ""
    curseqid = ""
    addlabel = False
    addbool = True
    onlyadd = False

    for char in treeline[::-1]:
        if onlyadd:
            newtreeline += char
            continue
        if char == "'":

            if not seqidbool:
                seqidbool = True
                newtreeline += char
                continue
            else:
                seqidbool = False
                seqidcount += 1

                curseqid = curseqid[::-1]

                if curseqid == collseqid[0]:
                    newtreeline += char
                    onlyadd = True
                    continue

                curseqid = ""
        if char == "&":

            newtreeline = newtreeline[:newtreeline.rfind(",")][::-1]

            add = f"&nodelabel={nodelabel},"

            newtreeline = add + newtreeline
            newtreeline = newtreeline[::-1]
            continue

        if seqidbool:
            curseqid += char
        newtreeline += char

    return newtreeline[::-1]


def count_labels(seqids, seqid2anot, subindexbool, speciesanot):

    allseqids = []
    speciesset = set()
    idbool = False
    for seqid in seqids:
        if subindexbool:
            seqid = ".".join(seqid.split(".")[:-1])
            allseqids.append(seqid)
            ids = seqid2anot.id902all.get(seqid)
        else:
            ids = seqid2anot.id902all.get(seqid)

            if not ids:
                subseqid = ".".join(seqid.split(".")[:-1])
                ids = seqid2anot.id902all.get(subseqid)
                if ids:
                    allseqids.append(subseqid)
                    subindexbool = True
                else:
                    allseqids.append(seqid)
        if ids:
            idbool = True
            allseqids.extend(ids)

    query = False
    speciesids = []
    anots = []
    for seqid in allseqids:

        anot = seqid2anot.seqid2anot.get(seqid)

        if not anot:
            anot = seqid2anot.seqid2anot.get(".".join(seqid.split(".")[:-1]))
            if not anot:
                query = True
                continue
        prevlen = len(speciesset)
        species = anot[0][speciesanot]

        if species.startswith("Unclassified") or species.startswith("-"):
            speciesids.append(seqid)
            anots.append(anot[0])
            continue

        speciesset.add(anot[0][speciesanot])

        if prevlen == len(speciesset):
            continue
        else:
            speciesids.append(seqid)
            anots.append(anot[0])

    if len(speciesids) == 0:
        print("0000")
        print(allseqids)
        return False, False, False, False, False

    return recursive_label_count(speciesids, anots) + (subindexbool,)

def format_taxa_counts(thedict):
    outlist = []

    for k, v in thedict.items():

        outlist.append(f"{k}({v})")

    return ";".join(outlist)


def recursive_label_count(speciesids, anots):

    ncbilist = ["Ncbi_order", "Ncbi_class", "Ncbi_phylum"]
    gtdblist = ["Gtdb_order", "Gtdb_class", "Gtdb_phylum"]
    nfirstcolor, gfirstcolor = None, None
    gtaxa2count = {}
    ntaxa2count = {}
    if not curtaxonomy or curtaxonomy == "Ncbi":
        for nind, el in enumerate(ncbilist):
            nbool, ntaxa2count, nfirstcolor, ntaxa, ncount, mix = get_taxa_labels(el, anots, speciesids)
            if nbool:
                break
    if not curtaxonomy or curtaxonomy == "Gtdb":
        for ind, el in enumerate(gtdblist):
            if ind < len(gtdblist) -1:
                upper = gtdblist[ind + 1]
            else:
                upper = False

            gbool, gtaxa2count, gfirstcolor, gtaxa, gcount, gmix = get_taxa_labels(el, anots, speciesids, upper=upper)

            if gbool:
                break

    if curtaxonomy == "Ncbi":
        firstcolor = nfirstcolor
    else:
        firstcolor = gfirstcolor

    ind, nind = 3,3
    ntaxa2count = format_taxa_counts(ntaxa2count)
    gtaxa2count = format_taxa_counts(gtaxa2count)

    if not ntaxa2count:
        label = f"gtdb_{gtaxa2count}"
        return {}, gtaxa2count, firstcolor, label
    elif not gtaxa2count:
        label = f"ncbi_{ntaxa2count}"

        return ntaxa2count, {}, firstcolor, label
    if ind < nind and gtaxa and  not gtaxa.startswith("Unclass"):
        label = f"gtdb_{gtaxa2count}"
        return ntaxa2count, gtaxa2count, gfirstcolor, label
    if gtaxa and gtaxa.startswith("Unclass"):
        label = f"ncbi_{ntaxa2count}"
    elif (ncount - gcount)/len(speciesids) > 0.3:

        label = f"ncbi_{ntaxa2count}"

    elif gmix and not mix:
        label = f"ncbi_{ntaxa2count}"

    else:
        label = f"gtdb_{gtaxa2count}"

    if not label:
        raise Exception("No label derived for collapsed clade")

    return ntaxa2count, gtaxa2count, firstcolor, label

def get_taxa_labels(el, anots, speciesids, **kwargs):

    taxa2count = dict()
    totalcount = 0
    taxa2anots = dict()

    iupper = kwargs.get("upper", el)
    if not iupper:
        iupper = el
    for ind, species in enumerate(speciesids):

        anot = anots[ind]
        taxa = anot[el]

        taxa2anots[taxa] = [anot["Colormap"], anot[iupper]]
        taxa2count[taxa] = taxa2count.get(taxa, 0) + 1
        totalcount += 1
    taxa2count = {k: v for k, v in sorted(taxa2count.items(), key=lambda item: item[1], reverse=True)}
    firstcounts = 0
    firstcolor = "#D3D3D3"
    firsttaxa = ""
    nbool = False
    newtaxa2count = dict()
    for taxa, counts in taxa2count.items():

        if len(taxa) == 1:
            taxa = "Unclassified"
        if taxa[1].isupper() or taxa[1].isdigit():

            curanot = taxa2anots[taxa]
            uppertaxa = curanot[1]

            taxa = f"{uppertaxa}_{taxa}"
            taxa2anots[taxa] = curanot
        if counts/totalcount >= 0.7 and not taxa.startswith("Unclass"):
            nbool = True
            newtaxa2count[taxa] = f"{counts}/{totalcount}"

        elif taxa.startswith("Unclass") and firstcounts:

            if (counts + firstcounts)/totalcount >= 0.8:
                nbool = True
                if counts/totalcount > 0.1 and counts > 1:
                    newtaxa2count[taxa] = f"{counts}/{totalcount}"
        elif counts/totalcount > 0.1 and counts > 1:

            newtaxa2count[taxa] = f"{counts}/{totalcount}"

        elif not firstcounts:
            newtaxa2count[taxa] = f"{counts}/{totalcount}"

        if not firstcounts:
            firstcounts = counts
            firsttaxa = taxa
            try:
                firstcolor = taxa2anots[taxa][0]
            except KeyError:
                firstcolor = "#808080"
    if firstcolor == "#000000":
        firstcolor = "#808080"

    if firstcounts/totalcount < 0.7:
        mix = f"({totalcount - firstcounts}/{totalcount})"
    else:
        mix = False

    return nbool, newtaxa2count, firstcolor, firsttaxa, firstcounts, mix

def  parse_nexus_to_collapse(treefile, seqidanot, attrlines, reanotbool, taxonomy):



    with open(treefile) as tree:

        format = next(tree)

        if format.strip() != "#NEXUS":
            raise Exception(f"The given file {treefile} is not in nexus format. "
                            f"Check your input directory and extension given for reanotation")
        outpath = f"{os.path.splitext(treefile)[0]}_collapse.figtree"
        outfile = open(outpath, "w")
        outfile.write(format)
        start = False
        anotbool = True
        itree = 1
        for line in tree:

            if line.startswith("begin trees"):
                start = False
                outfile.write(line)
            elif line.strip().startswith(";") and start:
                start = False
                outfile.write(line)

            elif start and reanotbool:
                seqid = line.split("[")[0]
                seqidbase = seqid.strip().replace("'", "")

                curseqid, lineadd = seqidanot.func(seqidbase, anotbool)
                if not curseqid.endswith("\n"):
                    curseqid = f"{curseqid}\n"
                outfile.write(curseqid)

            elif start and not reanotbool:
                outfile.write(line)
            elif line.strip().startswith("tree tree"):

                newtreeline = treeline_collapse(line, seqidanot, itree, os.path.splitext(treefile)[0], taxonomy)
                outfile.write(newtreeline[1])
            elif line.strip().startswith("taxlabels"):
                start = True
                outfile.write(line)
            elif line.strip().startswith("begin figtree"):
                break
            else:
                outfile.write(line)


        for line in attrlines:
            outfile.write(line)

        return outpath

def reanotate_nexus(seqidanot, treefile, reanotate, taxonomy):

    sys.stdout.write(f"Start collapsing {treefile}\n")
    taxlines, attrlines = parse_figtree_config("nexusfigtree_set_to_collapse.txt")
    return parse_nexus_to_collapse(treefile, seqidanot, attrlines, reanotate, taxonomy)




