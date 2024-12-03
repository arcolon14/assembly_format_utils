# Assembly format utilities

Different utility scripts for formatting and organizing genome assemblies and annotations.

## Rename and sort genome assemblies

```sh
$ python3 rename_sort_fa.py -h
  usage: rename_sort_fa.py [-h] -f IN_FASTA -i IN_FAI [-k NAME_KEY] [-o OUT_DIR] [-b BASENAME] [-m MIN_CHR_LEN] [-g GTF] [-l MIN_SEQ_LEN]
                           [--rename-by-length] [--export-sorted]
  
  Process an input FASTA along with annotations to rename the sequences and sort the output.
  
  options:
    -h, --help            show this help message and exit
    -f IN_FASTA, --in-fasta IN_FASTA
                          (str) Path to input FASTA.
    -i IN_FAI, --in-fai IN_FAI
                          (str) Path to input FASTA index (FAI).
    -k NAME_KEY, --name-key NAME_KEY
                          (str) Path to name key file.
    -o OUT_DIR, --out-dir OUT_DIR
                          (str) Output directory.
    -b BASENAME, --basename BASENAME
                          (str) Name of current run. Defaults to datetime.
    -m MIN_CHR_LEN, --min-chr-len MIN_CHR_LEN
                          (int/float) Minimum length of sequence to label as a chromosome. [default = 1,000,000]
    -g GTF, --gtf GTF     (str) Path to input GTF/GFF3.
    -l MIN_SEQ_LEN, --min-seq-len MIN_SEQ_LEN
                          (int/float) minimum length needed to export a sequence. [default = 1,000]
    --rename-by-length    Rename the chromosome sequences by their length in BP (i.e., longest sequences is chromosome 1).
    --export-sorted       Export the sequences sorted by size. Defaults to the order of sequences in the FAI.
```

### Rename key example

```sh
# OldName    NewName
scaffold_2   01
scaffold_6   02
scaffold_5   03
scaffold_7   04
scaffold_8   05
scaffold_9   06
```

## Author

Angel G. Rivera-Colon  
Institute of Ecology and Evolution  
University of Oregon

