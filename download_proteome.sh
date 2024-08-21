#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH --job-name=download #job name
#SBATCH -D /projects/tollis_lab/mammals_gene_dups/mammals_proteomes # set working directory to

module load anaconda3/2023.07
conda activate mammals

#conda install -c bioconda ncbi-genome-download

#species with available annotations listed in available_annotations.txt
#get list of accessions
cut -f 2 available_annotation.txt > available_accessions.txt

#download annotations from GCA accessions
while read accession; do
ncbi-genome-download -A "$accession" -s genbank -F protein-fasta --flat-output -o ./ vertebrate_mammalian
done < /projects/tollis_lab/mammals_gene_dups/mammals_annotations/available_accessions.txt

#download annotaions from GCF accessions (refseq)
while read accession; do
ncbi-genome-download -A "$accession" -F protein-fasta --flat-output -o ./ vertebrate_mammalian
done < /projects/tollis_lab/mammals_gene_dups/mammals_annotations/available_accessions.txt
