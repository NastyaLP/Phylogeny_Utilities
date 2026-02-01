import sys, os

def get_round(label):

    ufb_color = {90: "#000000", 75: "#414B56", 50: "#A8ADB4", 0: "#FFFFFF"}

    try:
        bootstrap = int(label)
    except TypeError:
        print(f"Not a bootstrap int to round {label}")
        raise
    if bootstrap < 50:
        return 0,  ufb_color[0]
    elif 50 <= bootstrap < 75:
        return 50, ufb_color[50]
    elif 75 <= bootstrap < 90:
        return 75, ufb_color[75]
    else:
        return 90, ufb_color[90]

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


                newtreeline += curseqid
                seqids.append(curseqid)

            else:
                if labelbool:
                    labelbool = False
                    label_round, color_round = get_round(newtreeline[-3:].lstrip("l").lstrip("="))
                    newtreeline += f", ufb_rounded={label_round}, ufb_color={color_round}"
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

    #print(newtreeline)
    return seqids,  newtreeline
