# bioinfo-scripts
Simple, this repository contains useful scripts to manipulate and visualize genomic data.


[TOC]: #

# Table of Contents
- [Convert a variant (SNP) table from genotypes to binary table](#convert-variant-table-from-genotypes-to-binary-table)



## Convert a variant (SNP) table from genotypes to binary table
[`snp_genotypes_to_binary_table.py`](snp_genotypes_to_binary_table.py)

```
Usage: snp_genotypes_to_binary_table.py [-h] -i INPUT -o OUTPUT

Converts a comma-separated SNP table as the output of RedDog pipeline, from genotypes (nucleotide) data into a binary table, taking into account a "Reference" column
necessarily present in the input table. Genotypes equal to the "Reference" column or gaps ("-") will be attributed a "0" value, whereas genotypes that differ from the
"Reference" one will be attributed a "1" value.

Optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        name of the input comma-separated genotypes SNP table
  -o OUTPUT, --output OUTPUT
                        name of the output comma-separated binary SNP table

						 
```

## Calculate sequencing depth by windows given by a BED file
[`calculate_averages_depths_bed_windows.py`](calculate_averages_depths_bed_windows.py)

```
Usage: calculate_averages_depths_bed_windows.py [-h] -i DEPTHS -b BED -o OUTPUT

Calculates the average sequencing depths in windows from a comma-separated table like the output of bedtools genomecov or multicov, given an imput BED file with the start and end positions of each window.

  -h, --help            show this help message and exit
  -i DEPTHS, --depths DEPTHS
                        Filename of the table of depths like given by bedtools, in
                        comma-separated format. First column correspond to the reference
                        chromosome position, second column corresponds to the sequencing
                        depth
  -b BED, --bed BED     Filename of the BED file. The BED file is a tab-delimited file
                        with three columns: a first column details the name of the
                        reference chromosome, second and third columns are the start and
                        end of each window
  -o OUTPUT, --output OUTPUT
                        Filename of the output tab-delimited file

```
