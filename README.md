**Anotate phylogenies**

To see options:

python3 convert_anotate_to_pdf_domain.py --help



**Usage Examples (examples for reproducing commands provided within test_dir)**


1. **To anotate and color and round bootstraps for plain newick treefile:**

python3 convert_anotate_to_pdf_domain.py test_dir rooted --annotation test_dir/anotfile.tsv --convert

where:

- "test_dir" for the directory with trees of any depth 
- "rooted" for the file extensions to look  

! In anotfile required columns are Seqid (for sequence id) and Colormap (for colors you want to map according to selected anotation feature)


2. **To collapse cartoon**

python3 convert_anotate_to_pdf_domain.py test_dir cartoon.figtree --collapse --annotation test_dir/anotfile.tsv --reanotate --convert

where:

"cartoon.figtree" Extensions of your cartooned phylogenies that you want to process for clade labels and colors
"reanotate" Reanotate tree (in case there were changes in annotation file)

You need additional columns for that option in your anotation file:
    I will write if the need exists