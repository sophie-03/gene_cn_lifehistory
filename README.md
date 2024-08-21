# The role of gene copy number in adaptation to life history traits across mammals 

## Repository information

This repository contains the scripts for analysis carried out in the paper "The role of gene copy number in adaptation to life history traits across mammals".

## Workflow

1. Data Collection
2. Identification of orthogroups
3. PGLS (individual genes)
4. Aggregate PGLS (gene sets)
5. Randomisation test
9. GSEA
10. ORA

### 1. Data Collection
This study uses 105 species, of these 105 species, 54 have available proteomes. These proteomes are downloaded from NCBI using [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download), in script ```download_proteome.sh```.

For the remaining species, proteomes were created from annotations using [gfftk](https://github.com/nextgenusfs/gfftk/tree/main) (v23.11.2), in script ```gfftk_proteome.sh```.

Keep only the longest transcript for each protein, using the [Orthofinder](https://github.com/davidemms/OrthoFinder) script [primary_transcript.py](https://github.com/davidemms/OrthoFinder/blob/master/tools/primary_transcript.py).

BUSCO is run to assess the quality of proteomes, in script ```busco_proteome.sh```.

### 2. Identification of Orthogroups
Gene orthogroups were inferred using [Orthofinder](https://github.com/davidemms/OrthoFinder) (v2.5.5)​, to estimate the gene copy number for all protein coding genes, in script ```orthofinder.sh```.

Species with low BUSCO scores were removed from analysis in ```orthofinder_remove_species.sh```.

Identification of which genes corresponded to each orthogroup was done using house mouse (*Mus musculus*) gene annotations, in script ```mouse_orthogroup_genes.sh```.

### 3. PGLS (individual genes)
PGLS to test for association between gene copy number and life history traits. Scripts in ```pgls_indv``` folder.

### 4. Aggregate PGLS (gene sets)
PGLS to test for association between the aggregate copy number of genes within gene sets and life history traits. Scripts in ```pgls_sets``` folder.

### 5. Randomisation test
Randomisation test: given a set of interest, produce 1000 replicates of the set by randomly selecting genes with a similar variance to those found within the original set. Carry out PGLS on the replicated sets to compare test statistics to original results, in script ```pgls_sets_simulation.R```.
