#!/usr/bin/env bash
#

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <dictionary_file> <csv_file> <output_file>"
    exit 1
fi

# Command-line arguments
dictionary_file="$1"
csv_file="$2"
output_file="$3"

# Backup file
backup_file="${csv_file%.csv}_backup.csv"

# Create a backup of the original CSV file
cp "$csv_file" "$backup_file"

# Remove the carriage return characters inside the dictionary file
sed -i 's/\r//' "$dictionary_file"

# Iterate through the lines of the dictionary file
while IFS=$'\t' read -r old_sample_id new_sample_id; do
    # Use awk to replace old sample IDs with new ones in the CSV file
    awk -v old="$old_sample_id" -v new="$new_sample_id" 'BEGIN {FS=OFS=","} {for (i=1; i<=NF; i++) if ($i == old) $i = new} 1' "$csv_file" > temp.csv
    mv temp.csv "$csv_file"
done < "$dictionary_file"

echo "Replacement completed. Output file: $output_file"
echo "Original data backed up to: $backup_file"