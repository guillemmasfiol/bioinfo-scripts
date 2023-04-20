#!/usr/bin/env python3
# **** Description***
# Calculates average sequencing depths from a comma-separated table like the output of bedtools genomecov or multicov, given an imput BED file with the start and end positions of each window.
#
# Author(s)
# Guillem Mas Fiol (guillem.mas-fiol@pasteur.fr) 



import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='Calculate average depths in windows of a reference chromosome given by a BED file. Both the input table and BED file cannot contain column names.')
parser.add_argument('-i', '--depths', required=True,  help='Filename of the table of depths like given by bedtools, in comma-separated format. First column correspond to the reference chromosome position, second column corresponds to the sequencing depth')
parser.add_argument('-b', '--bed', required=True,  help='Filename of the BED file. The BED file is a tab-delimited file with three columns: a first column details the name of the reference chromosome, second and third columns are the start and end of each window')
parser.add_argument('-o', '--output', required=True,  help='Filename of the output tab-delimited file')
args = parser.parse_args()


# Read input
depths_file = args.depths
depths = pd.read_csv(depths_file, sep=',', header=None, names=['Position', 'Depth'])

# Read the BED file defining the windows
bed_file = args.bed
bed = pd.read_csv(bed_file, sep='\t', header=None, names=['Contig', 'Start', 'End'])


# Function to calculate the average depth in each window
def calculate_window_depth(row, depths):
    start, end = float(row['Start']), float(row['End'])
    depths_in_window = depths[(depths['Position'] >= start) & (depths['Position'] < end)]
    return pd.Series({
        'Contig': row['Contig'],
        'Start': start,
        'End': end,
        'Depth': depths_in_window['Depth'].mean()
    })

# Apply the function to each row of the BED file and store the results in a new DataFrame
window_depths = bed.apply(calculate_window_depth, axis=1, depths=depths).reset_index(drop=True)


# Save output to file
window_depths.to_csv(args.output, sep='\t', index=False)

