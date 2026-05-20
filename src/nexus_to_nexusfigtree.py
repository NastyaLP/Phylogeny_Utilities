import os, sys
import logging



class SeqidAnot:

    def __init__(self, anotfile, colordict, column):

        self.func = None
        self.seqid2anot = dict()
        self.enzyme2count = dict()
        self.anot2count = dict()
        self.seqid2subindex = dict()
        self.id902all = dict()
        self.enzbool = True
        self.parse_anot(anotfile, colordict, column)

        self.merge = False

    def set_merge(self, merge, merge_annotations):

        self.merge = merge

        if not self.merge:
            return False

        if merge_annotations:

            merge_annotations = merge_annotations.split(",")

            for el in merge_annotations:

                self.anot2count[el] = dict()

        else:

            anotdict = next(iter(self.seqid2anot.values()))

            for key, value in anotdict[0].items():

                self.anot2count[key] = dict()

    def reset_merge(self):

        self.enzyme2count = dict()

        for key, value in self.anot2count.items():

            self.anot2count[key] = dict()

    def get_anot_func(self, seqid, anotbool):

        anotdict = self.seqid2anot.get(seqid)

        if not anotdict:

            try:
                newseqid = ".".join(seqid.split(".")[:-1])
                anotdict = self.seqid2anot[newseqid]
                curid = self.seqid2subindex.get(newseqid, 0)

                if anotbool:
                    curid += 1
                    seqid = f"{newseqid}.{curid}"
                    self.seqid2subindex[newseqid] = curid
                seqid = f"{newseqid}.{curid}"

            except KeyError:
                if "." in seqid:
                    newseqid = ".".join(seqid.split(".")[:-1])
                    curid = self.seqid2subindex.get(newseqid, 0)
                    if anotbool:
                        curid += 1
                        self.seqid2subindex[newseqid] = curid
                    seqid = f"{newseqid}.{curid}"
                logging.warning(f"Seqid '{seqid}' is not in anotfile will be colored in black")
                lineadd = self.construct_treeline("#d3d3d3", {})
                return f"'{seqid}'", lineadd
        seqid = self.construct_seqid(seqid, anotdict[0])
        lineadd = self.construct_treeline(anotdict[1], anotdict[2])

        if self.enzbool:

            self.update_enz_counts(anotdict[2]["enzyme"])

        if self.merge:

            self.update_counts(anotdict[0])

        return seqid, lineadd

    def update_enz_counts(self, enzyme):

        count = self.enzyme2count.get(enzyme, None)

        if not count:
            self.enzyme2count[enzyme] = 1

        else:
            count += 1
            self.enzyme2count[enzyme] = count

    def update_counts(self, anotdict):

        for key, value in self.anot2count.items():

            curanot = anotdict[key]

            count = value.get(curanot)
            if not count:
                value[curanot] = 1

            else:
                value[curanot] += 1

    def construct_seqid(self, seqid, anotdict):

        anotline = f"'{seqid}'[&"
        anotlinelist = []
        for key, value in anotdict.items():

            value = f'"{value}"'
            anotlinelist.append("=".join([key, value]))

        anotline += ",".join(anotlinelist)
        anotline += "]\n"

        return anotline

    def construct_treeline(self, color, addanotdict):

        treeline = f"'[&!color={color}, thick=60"

        for key, value in addanotdict.items():

            treeline += f", {key}={value}"

        return treeline + "]"

    def get_seqid_func(self, seqid):

        return f"'{seqid}'", "'"

    # for color anot to work colormap column should be present, enzyme column for legend and shapes is optional
    def parse_anot(self, anotfile, colordict, column):
        id90bool = True
        if not anotfile:

            self.func = self.get_seqid_func
            return False
        with open(anotfile, "r") as inanot:
            head = next(inanot).strip("\n").split("\t")

            iseqid = head.index("Seqid")

            head.pop(iseqid)

            try:
                id90 = head.index("Id90")
            except ValueError:
                id90bool = False
                pass
            try:
                ienzyme = head.index("Enzyme")

                def lineextract(curline):
                    return {"enzyme": curline[ienzyme]}

            except ValueError:
                ienzyme = False
                self.enzbool = False

                def lineextract(curline):
                    return {}


            try:
                icolumn = head.index(column)
            except ValueError:
                try:
                    icolumn = head.index(column.lower())
                except ValueError:
                    logging.error(f"The column {column} is absent. Either provide correct name of the column with"
                                  f"--color_by_column or have default column 'Gtdb_domain' in your annotations")
                    exit()

            anotlist = [el for el in head]

            for line in inanot:

                line = line.strip("\n").split("\t")

                seqid = line[iseqid]
                line.pop(iseqid)

                color = colordict.get(line[icolumn])

                enzyme = lineextract(line)
                if id90bool:
                    if line[id90]:
                        self.id902all[line[id90]] = self.id902all.get(line[id90], []) + [seqid]

                anotline = {el: line[iel] for iel, el in enumerate(anotlist)}

                self.seqid2anot[seqid] = [anotline, color, enzyme]

        self.func = self.get_anot_func


# for color anot to work colormap column should be present
#simple nexus or newick as an input
def parse_nexus_newick(infile, seqidanot, splitbool):

    outre = None
    lines = []
    with open(infile) as inf:
        try:
            firstline = next(inf)
        except:
            logging.error(f"{infile} is not newick/nexus format")
            exit(1)

        if firstline.startswith("("):
            if splitbool:
                outre, lines = newick_parser_split(firstline, inf, seqidanot, infile)

            else:
                outre, lines = newick_parser(firstline,  inf, seqidanot, infile)

        elif firstline.startswith("#NEXUS"):
            #if splitbool:

            #    outre, lines = nexus_parser_split(firstline, inf, seqidanot) not implemented
            #else:

            outre, lines = nexus_parser(firstline,  inf, seqidanot, infile)

        else:
            raise Exception("Unknown tree format")

        return outre, lines


def newick_parser(firstline, inf, seqidanot, infile):

    n = 0
    insert = "tree tree{iin} = [&R] "
    anotbool = True
    lines = [f"Begin trees;\n"]
    outre, newline = get_id_nexus(firstline, seqidanot, infile, anotbool)
    newline = insert.format(iin=n) + newline

    lines.append(newline)
    n += 1
    anotbool = False
    for line in inf:
        if line.startswith("("):
            outre, newline = get_id_nexus(line, seqidanot, infile,anotbool)
            newline = insert.format(iin=n) + newline
            lines.append(newline)
            n += 1
    lines.append(f"end;\n")
    return outre, lines


# to the split files
def newick_parser_split(firstline, inf, seqidanot, infile):

    insert = "tree tree0 = [&R] "
    start = f"Begin trees;\n"
    end = f"end;\n"
    totalines = []

    outre, newline = get_id_nexus(firstline, seqidanot, infile)

    newline = insert + newline
    totalines.append([start, newline, end])

    for line in inf:

        if line.startswith("("):
            outre, newline = get_id_nexus(line, seqidanot, infile)
            newline = insert + newline
            totalines.append([start, newline, end])
    return outre, totalines


def nexus_parser(firstline, inf, seqidanot, infile):
    lines = []
    print(firstline)
    for line in inf:

        if line.startswith("tree"):
            outre, newline = get_id_nexus(line, seqidanot, infile)
            lines.append(newline)
        else:
            lines.append(line)

    return outre, lines


def parse_figtree_config(configname):

    figtreeconfig = os.path.join(os.path.dirname(__file__), configname)

    taxlines = []
    attrlines = []

    attrbool = False
    with open(figtreeconfig) as infigc:

        for line in infigc:

            if line.startswith("##"):
                attrbool = True
                continue

            if not attrbool:
                taxlines.append(line)
            else:
                attrlines.append(line)

    return taxlines, attrlines


def write_nexus_figtree(outre, lines, taxlines, attrlines, outfile):

    ntax = len(outre)

    for line in taxlines:

        if line.lstrip().startswith("dimensions"):
            line = line.split("=")[0]
            line = "=".join([line, str(ntax) + ";" + "\n"])
            outfile.write(line)
            continue
        elif line.lstrip().startswith("taxlabels"):
            outfile.write(line)
            for el in outre:
                outfile.write("\t" + el)
            continue
        outfile.write(line)

    outfile.writelines(lines)
    
    outfile.writelines(attrlines)

def get_round(label):

    ufb_color = {90: "#000000", 75: "#414B56", 50: "#A8ADB4", 0: "#FFFFFF"}

    try:
        bootstrap = int(label)
    except TypeError:
        print(f"Not a bootstrap int to round {label}")
        raise
    if bootstrap < 50:
        return 0,  ufb_color[0]
    elif 50 <= bootstrap < 90:
        return 50, ufb_color[50]
    else:
        return 90, ufb_color[90]

def compare_ids(curseqid, newtreeline):

    treesuffix = newtreeline[-6:]
    if "." not in treesuffix:
        return newtreeline

    treesuffix = treesuffix.split(".")[-1]
    cursuffix = curseqid.split("[")[0].strip("'").split(".")[-1]

    if cursuffix != treesuffix:

        newtreeline = newtreeline[:-len(treesuffix)] + cursuffix

    return newtreeline

def get_id_nexus(treeline, seqidanot, infile,anotbool):

    seqids = []
    newtreeline = ''
    openbool = False
    seqidbool = False
    closebool = False
    curseqid = ""

    seqidanot.reset_merge()
    countroot = 0
    seqidcount = 0
    countbool = False
    rootbool = False
    labelbool = False
    for char in treeline:

        if char == "(":
            if not countbool:
                countroot += 1
            openbool = True
            seqidbool = True
            newtreeline += char
            continue

        elif char == ")":

            if not countbool:
                countroot -=1
                if countroot == 1:
                    countbool = True
                    rootbool = True

            seqidbool = False
            openbool = False
            closebool = True
            newtreeline += char
            continue

        elif char == ":":
            if seqidbool:
                seqidbool = False
                seqidcount += 1

                curseqid, lineadd = seqidanot.func(curseqid, anotbool)

                newtreeline = compare_ids(curseqid, newtreeline)

                newtreeline += lineadd
                seqids.append(curseqid)

            else:
                if labelbool:
                    labelbool = False
                    label_round, color_round = get_round(newtreeline[-3:].lstrip("l").lstrip("="))
                    newtreeline += f", ufb_rounded={label_round}"
                else:
                    if rootbool:
                        newtreeline += "[&!color=#aa4a44, thick=100]"
                        rootbool = False
                    else:
                        newtreeline += "[&thick=60"

                newtreeline += "]"
                closebool = False
            newtreeline += char
            continue

        elif char == ",":

            newtreeline += char
            #newtreeline += "'"
            openbool = True
            curseqid = ""
            seqidbool = True
            continue

        if seqidbool:
            if openbool:
                newtreeline += "'"
                openbool = False
            curseqid += char
        elif closebool and not seqidbool and char != ";":
            closebool = False
            if rootbool:
                newtreeline += "[&!color=#aa4a44, thick=100, "

                rootbool = False

            else:
                newtreeline += "[&thick=60, "

            newtreeline += "label="
            labelbool = True
        newtreeline += char
        if char == ";":
            newtreeline += "\n"
            break
    if countbool and rootbool:
        print(f"Single brunch rooting {infile}, manual attention recommended")
        newtreeline = single_brunch_rooting(newtreeline)

    return seqids,  newtreeline

def single_brunch_rooting(newtreeline):

    rootedtreeline = ""
    rootbool = False
    for char in newtreeline:

        if char == ":" and not rootbool:

            templine = rootedtreeline.split(", ")
            newtempline = []
            for el in templine:

                if el.startswith("thick"):
                    el = "thick=100"
                newtempline.append(el)
            rootedtreeline = ", ".join(newtempline)
            rootbool = True

        rootedtreeline += char

    return rootedtreeline



