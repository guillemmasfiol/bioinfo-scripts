# bioinfo-scripts
Simple, this repository contains useful scripts to manipulate and visualize genomic data.


[TOC]: #

# Table of Contents
- [Plasmid2Chromosome Depth Ratios Calculator](#Plasmid2Chromosome-Depth-Ratios-Calculator)
- [Convert a variant (SNP) table from genotypes to binary table](#convert-variant-table-from-genotypes-to-binary-table)
- [Calculate sequencing depth by windows given by a BED file](#calculate-sequencing-depth-by-windows-given-by-a-BED-file)
- [Extracts annotations from genbank file within a given range](#extracts-annotations-from-genbank-file-within-a-given-range)

## Plasmid2Chromosome Depth Ratios Calculator
[`plasmids2chromosome_ratio_calculator.sh`](plasmids2chromosome_ratio_calculator.sh)
```
Calculates the plasmids:chromosome ratio using the depths of reads mapped to reference replicons for the chromosome and different given plasmid identifiers.
BAM files are stored in the same directory to automatically compute their sequencing depths.
It uses samtools depth command, so it requires previous installation of samtools dependencies

Usage: ./plasmids2chromosome_ratio_calculator.sh [-d <bam_dir>] [-o <output_file>] [-c <chr_replicon ID>] [-p1 <plasmid1_replicon ID>] [-p2 <plasmid2_replicon ID>] [-p3 <plasmid3_replicon ID>] [-h]

Necessary arguments:
-d <bam_dir>             Path to directory containing the indexed BAM files
-o <output_file>         Path to output file
-c <chr_replicon>        ID of the chromosome replicon used for mapping the reads
-p1 <plasmid1_replicon>  ID of the first plasmid replicon used for mapping the reads
-p2 <plasmid2_replicon>  ID of the second plasmid replicon used for mapping the reads
-p3 <plasmid3_replicon>  ID of the third plasmid replicon used for mapping the reads
 
```

## Convert a variant (SNP) table from genotypes to binary table
[`snp_genotypes_to_binary_table.py`](snp_genotypes_to_binary_table.py)

```
Usage: python snp_genotypes_to_binary_table.py [-h] -i INPUT -o OUTPUT

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
Usage: python calculate_averages_depths_bed_windows.py [-h] -i DEPTHS -b BED -o OUTPUT

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

## Extracts annotations from genbank file within a given range
[`extract_annotations_range_genbank.py`](extract_annotations_range_genbank.py)

```
Usage: python extract_annotations_range_genbank.py <genbank_file> <start_position> <end_position> <output_file>

Extracts the annotations for different fields, for elements within a range defined as input arguments (start and end positions) and export results into a table with each column corresponding to a different annotation field.

```
