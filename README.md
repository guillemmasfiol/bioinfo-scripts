# bioinfo-tools
Simple, this repository contains useful scripts to manipulate and visualize genomic data.


[TOC]: #

# Table of Contents
- [Convert a variant (SNP) table from genotypes to binary table](#convert-variant-table-from-genotypes-to-binary-table)



## Convert a variant (SNP) table from genotypes to binary table
[`python snp_genotypes_to_binary_table.py`][python snp_genotypes_to_binary_table.py]

```
Usage: snp_genotypes_to_binary_table.py [-h] -i INPUT -o OUTPUT

Converts a comma-separated SNP table as the output of RedDog pipeline, from genotypes (nucleotide) data into a binary table, taking into account a "Reference" column
necessarily present in the input table. Genotypes equal to the "Reference" column or gaps ("-") will be attributed a "0" value, whereas genotypes that differ from the
"Reference" one will be attributed a "1" value

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        name of the input comma-separated genotypes SNP table
  -o OUTPUT, --output OUTPUT
                        name of the output comma-separated binary SNP table

						 
```
