#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH --job-name=primary_transcripts #job name
#SBATCH -D /projects/tollis_lab/mammals_gene_dups/mammals_proteomes # set working directory to
#SBATCH --mem=150000
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.matthews5@universityofgalway.ie # email address
#SBATCH --out=primary_transcript.out
#SBATCH --time=2-00:00:00


# extract only the longest isoform (primary transcript) from the proteome

module load anaconda3/2023.07
conda activate placement

cd /projects/tollis_lab/mammals_gene_dups/mammals_proteomes/non_liftoff_species/proteomes

# use the primary_transcript.py script provided in orthofinder tools
for f in *faa ; do python ~/bin/orthofinder_tutorial/OrthoFinder/tools/primary_transcript.py $f ; done
