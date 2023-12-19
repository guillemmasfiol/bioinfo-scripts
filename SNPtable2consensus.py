#!/usr/bin/env python
#
# 
# **** Description***
# Generates consensus fasta sequences for each sample given a comma-separated SNP table and the corresponding reference genome fasta files
#
# Author(s) 
#	Guillem Mas Fiol (guillem.mas-fiol@pasteur.fr) 
#

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

def read_snp_table(file_path):
    """
    Read the SNP table and return a dictionary of positions and genotypes for each sample.
    """
    snp_data = {}
    with open(file_path, "r") as snp_file:
        header = snp_file.readline().strip().split(",")
        for line in snp_file:
            parts = line.strip().split(",")
            pos = int(parts[0])
            for i in range(2, len(parts)):
                sample = header[i]
                genotype = parts[i]
                if sample not in snp_data:
                    snp_data[sample] = {}
                snp_data[sample][pos] = genotype
    return snp_data

def build_consensus_sequence(reference_seq, snp_data):
    """
    Build a consensus sequence for each sample based on the SNP data and the reference sequence.
    """
    consensus_sequences = {}
    reference_length = len(reference_seq)

    for sample, snps in snp_data.items():
        consensus_sequence = []
        for pos in range(1, reference_length + 1):
            if pos in snps:
                consensus_sequence.append(snps[pos])
            else:
                consensus_sequence.append(str(reference_seq[pos - 1]))
        
        consensus_sequences[sample] = "".join(consensus_sequence)

    return consensus_sequences

def write_consensus_fasta(consensus_sequences, output_folder):
    """
    Write consensus sequences to separate FASTA files for each sample.
    """
    for sample, sequence in consensus_sequences.items():
        seq = Seq(sequence)
        record = SeqRecord(seq, id=sample, description="Consensus sequence for " + sample)
        output_file = f"{output_folder}/{sample}_consensus.fasta"
        with open(output_file, "w") as output_handle:
            SeqIO.write(record, output_handle, "fasta")

def main():
    parser = argparse.ArgumentParser(description="This script takes a comma-separated SNP table, typically the output from RedDog mapping and variant calling pipeline, and generates a consensus fasta sequence for each sample with the same length as the reference sequence, including the variants specific to each sample. Variants include ATGC genotypes and '-' for gaps. Reference fasta sequence must be the same used for calling variants in the SNP table.")

    parser.add_argument("--s", dest="snp_table", required=True, help="Input comma-separated SNP table")
    parser.add_argument("--r", dest="reference_fasta", required=True, help="FASTA reference sequence")
    parser.add_argument("--o", dest="output_folder", required=True, help="Output directory name")

    args = parser.parse_args()

    # Read reference sequence
    reference_records = SeqIO.to_dict(SeqIO.parse(args.reference_fasta, "fasta"))
    reference_id, reference_seq = list(reference_records.items())[0]

    # Read SNP table
    snp_data = read_snp_table(args.snp_table)

    # Build consensus sequences
    consensus_sequences = build_consensus_sequence(reference_seq.seq, snp_data)

    # Write consensus sequences to separate FASTA files
    write_consensus_fasta(consensus_sequences, args.output_folder)

if __name__ == "__main__":
    main()