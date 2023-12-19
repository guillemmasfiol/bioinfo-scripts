#!/usr/bin/env bash
#

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <dictionary_file> <tab_delimited_file> <output_file>"
    exit 1
fi

# Command-line arguments
dictionary_file="$1"
tab_delimited_file="$2"
output_file="$3"

# Backup file
backup_file="${tab_delimited_file%.tsv}_backup.tsv"

# Create a backup of the original tab-delimited file
cp "$tab_delimited_file" "$backup_file"

# Remove the carriage return characters inside the dictionary file
sed -i 's/\r//' "$dictionary_file"

# Iterate through the lines of the dictionary file
while IFS=$'\t' read -r old_sample_id new_sample_id; do
    # Use awk to replace old sample IDs with new ones in the tab-delimited file
    awk -v old="$old_sample_id" -v new="$new_sample_id" 'BEGIN {FS=OFS="\t"} {for (i=1; i<=NF; i++) if ($i == old) $i = new} 1' "$tab_delimited_file" > temp.tsv
    mv temp.tsv "$tab_delimited_file"
done < "$dictionary_file"

echo "Replacement completed. Output file: $output_file"
echo "Original data backed up to: $backup_file"
