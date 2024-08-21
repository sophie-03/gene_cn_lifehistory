#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH --job-name=gfftk #job name
#SBATCH -D /projects/tollis_lab/mammals_gene_dups/mammals_proteomes/proteomes/liftoff_species # set working directory to
#SBATCH --mem=150000
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.matthews5@universityofgalway.ie # email address
#SBATCH --out=gfftk.out
#SBATCH --time=2-00:00:00

module load anaconda3/2023.07
conda activate placement

# Read species names from liftoff_species.txt and process each species
while IFS= read -r species; do
    echo "Processing $species..."

    # Remove CDS features with length < 10 from the GFF file
    awk -F"\t" '$3 != "CDS" || (($3 == "CDS") && ($5 - $4 + 1) >= 10) || $1 ~ /^#/ {print}' \
    "/projects/tollis_lab/TE/mammals/genomes/5_BED_Annotations/$species/$species.gff" \
    > "intermediate_files/${species}_filtered.gff"

    # Create proteome using gfftk convert
    gfftk convert --input-format gff3 --output-format proteins \
    -i "intermediate_files/${species}_filtered.gff" \
    -f "/projects/tollis_lab/TE/mammals/genomes/3_ref_container/${species}.fa" \
    -o "${species}.faa"

    echo "Finished processing $species."
done < ../../liftoff_species.txt
