import pandas as pd
import argparse


parser = argparse.ArgumentParser(description='Convert a comma-separated SNP table as the output of RedDog pipeline, from genotypes (nucleotide) data into a binary table taking into account a "Reference" column necessarily present in the input table. Genotypes equal to the "Reference" column or gaps ("-") will be attributed a "0" value, whereas genotypes that differ from the "Reference" one will be attributed a "1" value')
parser.add_argument('-i', '--input', required=True, help='name of the input comma-separated genotypes SNP table')
parser.add_argument('-o', '--output', required=True, help='name of the output binary SNP table')
args = parser.parse_args()

input_file = args.input
mutation_table = pd.read_csv(input_file)

# Copy of the mutation table with only the 'Gene' and sample genotype columns
genotype_table = mutation_table[['Gene', 'Reference'] + list(mutation_table.columns[2:])].copy()

# Function to transform a genotype column into a binary column
def binarize_genotype(row):
    reference_genotype = row['Reference']
    sample_genotype = row[2:]
    return [1 if g != '-' and g != reference_genotype else 0 for g in sample_genotype]

# Apply the function to each row of the genotype table and expand the resulting lists into separate binary columns
genotype_table[list(genotype_table.columns[2:])] = genotype_table.apply(binarize_genotype, axis=1, result_type='expand')


# Save results
output_file = args.output
genotype_table.to_csv(output_file, index=False)