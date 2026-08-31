# Anotate phylogenies

## Installation

Python 3.8+  
Required Python packages:  
- biopython
- PyPDF2

1. Clone a repository:  

```bash
git clone https://github.com/NastyaLP/Phylogeny_Utilities  
```
2. Install dependencies with:

```bash
pip install -r Phylogeny_Utilities/requirements.txt
```
or manually:
```bash
pip install biopython PyPDF2
```
3. Download and copy figtree.jar to Phylogeny_Utilities

- download und unzip FigTree.v1.4.4.zip from GitHub with wget or manually  
```bash
wget https://github.com/rambaut/figtree/releases/download/v1.4.4/FigTree.v1.4.4.zip
unzip FigTree.v1.4.4.zip
```
- copy figtree.jar to Phylogeny_Utilities
```bash
cp FigTree\ v1.4.4/lib/figtree.jar Phylogeny_Utilities/
```

4. Switch to the  Phylogeny_Utilities and run phylogeny_util.py to see all available options:

```bash
python3 phylogeny_util.py --help
```
---

## Usage Examples 

**Examples for reproducing commands provided within directory test**


1. **Annotate, color, and round bootstrap values for plain Newick tree files**

```bash
python3 phylogeny_util.py "test" "rooted" --annotation test/anotfile.tsv --convert
```

where:

- "test" for the directory with trees of any depth 
- "rooted" for the file extensions to look  
- "--convert" for converting nexus to pdf image

! In anotfile required columns are Seqid (for sequence id) and Colormap (for colors you want to map according to selected anotation feature)

---

2. **Collapse cartooned clades and label them with taxonomic counts per clade**

Default:

```bash
python3 phylogeny_util.py "test" "cartoon.figtree" --collapse --convert
```
where:
- "test" for the directory with trees of any depth 
- "cartoon.figtree" Extensions of your cartooned phylogenies that you want to process for clade labels and colors
- "--convert" for converting nexus to pdf image

You need specific columns with taxonomy information in your annotation file to run collapse option on your cartooned phylogenies:  

- When using only Ncbi taxonomy annotation file must have four columns named:  
    - Ncbi_phylum, Ncbi_class, Ncbi_order, Ncbi_species

```bash
python3 phylogeny_util.py "test" "cartoon.figtree" --collapse --taxonomy Ncbi --convert
```

- When using only GTDB taxonomy annotation file must have four columns named:  
    - Gtdb_phylum, Gtdb_class, Gtdb_order, Gtdb_species


```bash
python3 phylogeny_util.py "test" "cartoon.figtree" --collapse --taxonomy Gtdb --convert
```

By default without "--taxonomy" option specified it will use both taxonomies. That means both the Ncbi and Gtdb columns mentioned above should be present in the annotation file.

---

3. **Advanced usage**

When sequences reduction was performed before phylogenetic reconstruction, then use the annotation file with representative id column 'Id90''(cluster representative)    
Having 'Id90' allows to count reduced sequences belonging to the clades of their representatives to fully represent taxonomic information.


```bash
python3 phylogeny_util.py "test" "cartoon.figtree" --collapse --annotation test/anotfile.tsv --convert
```