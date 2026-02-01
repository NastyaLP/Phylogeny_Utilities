import sys, os
import tempfile
from PyPDF2 import PdfReader, PdfWriter, PdfMerger
from PyPDF2.generic import AnnotationBuilder


def pdfmerge(anotdict, outmergedpath, anotbool, args):

    merger = PdfMerger()
    print("START PDF")

    for file_path in args:
        basefile = os.path.basename(os.path.splitext(file_path)[0])
        print(f"Merging 1st page from: {file_path}")

        # Only take the first page
        merger.append(file_path, pages=(0, 1))  # pages=(start, end) -> end is exclusive

        # Write out merged result
    with open(outmergedpath, "wb") as out:
        merger.write(out)


def pdfanotate(infile, anottext, writer, anotbool):

    with open(infile, 'rb') as inf:

        reader = PdfReader(inf)
        writer.add_page(reader.pages[0])


        if anotbool:
            curnumber = len(writer.pages) - 1

            # Create the annotation and add it
            annotation = AnnotationBuilder.text(rect=(45, 45, 45, 45),
                                            text=anottext, open=True, flags=0)
            writer.add_annotation(page_number=curnumber, annotation=annotation)


        return writer


def pdf_anotate_single(infile, anottext):

    writer = PdfWriter()

    pdfanotate(infile,anottext, writer, True)

    with open(infile, "wb") as output_stream:
        writer.write(output_stream)



