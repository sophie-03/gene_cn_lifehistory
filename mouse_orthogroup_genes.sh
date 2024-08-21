#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH --job-name=ortho_genes
#SBATCH -D /projects/tollis_lab/mammals_gene_dups/orthofinder/proteomes/
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.matthews5@universityofgalway.ie # email address
#SBATCH --time=2-00:00:00
#SBATCH --mem=150000

# Set the paths to your input files
orthogroups_file="/projects/tollis_lab/mammals_gene_dups/orthofinder/proteomes/mouse_orthogroups_accessions.txt"
gff_file="/projects/tollis_lab/TE/mammals/genomes/5_BED_Annotations/Mus_musculus/Mus_musculus.gff"

# Create a new file name
output_file="mouse_orthogroups_genes.txt"

# Process each line in the orthogroups file
while read -r orthogroup accession; do
    # Use grep to find entries in the gff file for the given accession number
    grep_result=$(grep -m 1 "Name=$accession" "$gff_file")

    # Extract the gene name from the grep result using awk
    gene_name=$(echo "$grep_result" | awk -F";" '{for(i=1;i<=NF;i++)if($i~/^gene=/)print $i}' | cut -d"=" -f2)

    # Print the original line with the new gene name to the output file
    echo "$orthogroup $accession $gene_name" >> "$output_file"
done < "$orthogroups_file"

echo "Processing complete. Results saved to $output_file"

#remove all extra lines
awk 'NF > 2' mouse_orthogroups_genes.txt > cleaned_mouse_orthogroups_genes.txt
mv cleaned_mouse_orthogroups_genes.txt mouse_orthogroups_genes.txt
