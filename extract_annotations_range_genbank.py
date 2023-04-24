#!/usr/bin/env python3
# **** Description***
# Extracts the annotations for different fields, for elements within a range defined as input arguments (start and end positions) and export results in a table with each column corresponding to a different annotation field.
#
# Author(s)
# Guillem Mas Fiol (guillem.mas-fiol@pasteur.fr) 



import sys
from Bio import SeqIO
import argparse


if len(sys.argv) != 5:
    print("Usage: python extract_annotations_range_genbank.py <genbank_file> <start_position> <end_position> <output_file>")
    sys.exit(1)


input_file = sys.argv[1]
start_position = int(sys.argv[2])
end_position = int(sys.argv[3])
output_file = sys.argv[4]

# Define the headers for each annotation field we want to extract
headers = ["start", "end", "type", "gene", "product", "locus_tag", "old_locus_tag", "protein_id", "EC_number", "inference", "GO_function", "GO_process"]


with open(output_file, 'w') as f:
    f.write("\t".join(headers) + "\n")
    

    for record in SeqIO.parse(input_file, "genbank"):
        for feature in record.features:
            if start_position <= feature.location.start.position <= end_position or start_position <= feature.location.end.position <= end_position:
            
                attributes = {"start": feature.location.start.position, "end": feature.location.end.position, "type": feature.type}
                for qualifier in feature.qualifiers:
                    if qualifier in headers:
                        attributes[qualifier] = feature.qualifiers[qualifier][0]
                
                f.write("\t".join(str(attributes.get(header, "")) for header in headers) + "\n")
