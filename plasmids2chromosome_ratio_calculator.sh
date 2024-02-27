#!/usr/bin/env bash

#
module load samtools
module load bwa
module load bedtools


function usage {
  echo ""
  echo ""
  echo "******************************************************************************"
  echo "***********     PLASMID2CHROMOSOME DEPTHS RATIO CALCULATOR       *************"
  echo "******************************************************************************"
  echo "*********	     Freely distributed 		    ******************"
  echo "******	      Share, cooperate, improve and transfer    **********************"
  echo "****** For the goodness of science, not for money and competition ************"
  echo "******************************************************************************"
  echo ""
  echo "Calculate depth ratios for specified replicons in a set of BAM files."
  echo ""
  echo "Author(s): Guillem Mas Fiol https://github.com/guillemmasfiol  guillem.mas-fiol@pasteur.fr"
  echo ""
  echo ""
  echo "Usage: $0 [-i <bam_dir>] [-o <output_file>] [-c <chr_replicon>] [-p1 <plasmid1_replicon>] [-p2 <plasmid2_replicon>] [-p3 <plasmid3_replicon>] [-h]"
  echo ""
  echo "Example: $0 -i bams_all -o depths.txt -c NC_003143.1 -p1 NC_003132.1 -p2 NC_003131.1 -p3 NC_003134.1"
  echo ""
  echo ""
  echo "Options:"
  echo "  -i <bam_dir>             Path to directory containing the indexed BAM files"
  echo "  -o <output_file>         Path to output file"
  echo "  -c <chr_replicon>        ID of the chromosome replicon used for mapping the reads"
  echo "  -p1 <plasmid1_replicon>  ID of the first plasmid replicon used for mapping the reads"
  echo "  -p2 <plasmid2_replicon>  ID of the second plasmid replicon used for mapping the reads"
  echo "  -p3 <plasmid3_replicon>  ID of the third plasmid replicon used for mapping the reads"
  echo "  -h                       Display this help message and exit"
  echo ""
  exit 1
}

while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -h|--help)
    usage
    exit 0
    ;;
    -i|--input_directory)
    bam_dir="$2"
    shift
    shift
    ;;
    -o|--output_file)
    output_file="$2"
    shift
    shift
    ;;
    -c|--chr_replicon)
    chr_replicon="$2"
    shift
    shift
    ;;
    -p1|--plasmid1_replicon)
    plasmid1_replicon="$2"
    shift
    shift
    ;;
    -p2|--plasmid2_replicon)
    plasmid2_replicon="$2"
    shift
    shift
    ;;
    -p3|--plasmid3_replicon)
    plasmid3_replicon="$2"
    shift
    shift
    ;;
	-h|--help)
      echo "$help_message"
      exit;;
    *)
    echo "Error: Invalid argument. Use --help to see the valid options"
    usage
    exit 1
    ;;
esac
done

# Check if all required input arguments are provided
if [ -z "$bam_dir" ] || [ -z "$output_file" ] || [ -z "$chr_replicon" ] || [ -z "$plasmid1_replicon" ] || [ -z "$plasmid2_replicon" ] || [ -z "$plasmid3_replicon" ]; then
  echo "Error: Missing required input arguments"
  usage
  exit 1
fi


bam_files=$(ls $bam_dir/*.bam)

# Calculate average read depth for each replicon across all samples
printf "Replicon\tAverage_Depth\n" > $output_file

for replicon in $chr_replicon $plasmid1_replicon $plasmid2_replicon $plasmid3_replicon; do
  depth=$(samtools depth -a -r $replicon $bam_files | awk '{sum+=$3} END {print sum/NR}')
  printf "%s\t%.2f\n" "$replicon" "$depth" >> $output_file
done

# Calculate the sequencing depth plasmid-to-chromosome ratios for each sample
printf "Sample\tChromosome_depth\tPlasmid1_ratio\tPlasmid2_ratio\tPlasmid3_ratio\n" > $output_file

for bam_file in $bam_files; do
  sample_id=$(basename $bam_file .bam)
  chr_depth=$(samtools depth -a -r $chr_replicon $bam_file | awk '{sum+=$3} END {print sum/NR}')
  plasmid1_depth=$(samtools depth -a -r $plasmid1_replicon $bam_file | awk '{sum+=$3} END {print sum/NR}')
  plasmid2_depth=$(samtools depth -a -r $plasmid2_replicon $bam_file | awk '{sum+=$3} END {print sum/NR}')
  plasmid3_depth=$(samtools depth -a -r $plasmid3_replicon $bam_file | awk '{sum+=$3} END {print sum/NR}')
  plasmid1_ratio=$(echo "$plasmid1_depth / $chr_depth" | bc -l)
  plasmid2_ratio=$(echo "$plasmid2_depth / $chr_depth" | bc -l)
  plasmid3_ratio=$(echo "$plasmid3_depth / $chr_depth" | bc -l)
  printf "%s\t%.2f\t%.2f\t%.2f\t%.2f\n" "$sample_id" "$chr_depth" "$plasmid1_ratio" "$plasmid2_ratio" "$plasmid3_ratio" >> $output_file
done
