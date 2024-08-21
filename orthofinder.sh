#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH --job-name=orthofinder1 #job name
#SBATCH -D /projects/tollis_lab/mammals_gene_dups/orthofinder # set working directory to
#SBATCH --mail-type=ALL # send email at job completion
#SBATCH --mail-user=s.matthews5@universityofgalway.ie # email address
#SBATCH --out=orthofinder.out
#SBATCH --time=7-00:00:00
#SBATCH --mem=150000
#SBATCH --cpus-per-task=32

module load anaconda3/2023.07
conda activate placement

## TUTORIAL
#cd example
#orthofinder
#orthofinder -f ExampleData/

cd /projects/tollis_lab/mammals_gene_dups/mammals_proteomes/proteomes/non_liftoff_species/primary_transcripts


#orthofinder
orthofinder -t 32 -f ./ -o /projects/tollis_lab/mammals_gene_dups/orthofinder/proteomes
