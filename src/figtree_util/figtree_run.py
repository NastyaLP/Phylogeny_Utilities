import os, sys



def figrun():

    figpath = f"{os.path.dirname(os.path.abspath(__file__))}/figtree.jar"

    teststr = ['java', '-jar', figpath, '-graphic', 'PDF', '-width', '800', '-height', '1020',
                                '/Users/stasia/PycharmProjects/treematrix/TreeVisual/figtree_util/test/test1_mfnF_auto.figtree', '/Users/stasia/PycharmProjects/treematrix/TreeVisual/figtree_util/test/test1_mfnF_auto.pdf']

    teststr = " ".join(teststr)
    print(teststr)
    os.system(teststr)
    sys.exit(1)
