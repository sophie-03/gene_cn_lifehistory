# The role of gene copy number in adaptation to life history traits across mammals 

## Repository information

This repository contains the scripts for analysis carried out in the paper "The role of gene copy number in adaptation to life history traits across mammals".

## Workflow

1. Creation of proteomes
2. BUSCO
3. Orthofinder
4. Labelling orthogroups with mouse genes
5. PGLS (individual genes)
6. Aggregate PGLS (gene sets)
7. Reconstructing TGFB
8. Simulation
9. GSEA
10. ORA

### 1. Creation of proteomes
This study uses 105 species, of these 105 species, 54 have available proteomes. These proteomes are downloaded from NCBI using [ncbi-genome-download](https://github.com/kblin/ncbi-genome-download), in script ```download_proteome.sh```.

For the remaining species, proteomes were created from annotations using [gfftk](https://github.com/nextgenusfs/gfftk/tree/main) (v23.11.2), in script ```gfftk_proteome.sh```.

Keep only the longest transcript for each protein, using the [Orthofinder](https://github.com/davidemms/OrthoFinder) script ```primary_transcript.py```.

### 2. BUSCO
BUSCO is run to assess the quality of proteomes, in script ```busco_proteome.sh```.
