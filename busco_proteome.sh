#!/bin/bash
#SBATCH --export=ALL # export all environment variables to the batch job
#SBATCH --job-name=busco #job name
#SBATCH --mem=150000
#SBATCH -D /projects/tollis_lab/mammals_gene_dups/mammals_proteomes/busco # set working directory to
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=s.matthews5@universityofgalway.ie # email address
#SBATCH --out=busco.out
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=6
#SBATCH --array=1-51

module load anaconda3/2023.07
conda activate placement

module load busco/5.4.7

# Define the input directory and output directory
input_dir="/projects/tollis_lab/mammals_gene_dups/mammals_proteomes/proteomes/liftoff_species/"

# Get the list of input files
files=(${input_dir}/*.faa)

# Get the current input file from the array index
input_file=${files[$SLURM_ARRAY_TASK_ID]}
filename=$(basename "$input_file")
filename_without_ext="${filename%.*}"

# Create a directory for each input file within the output directory
output_subdir="${filename_without_ext}"
mkdir -p "$output_subdir"

# Run the BUSCO command
busco -m protein -i "$input_file" -o "$output_subdir" -l mammalia_odb10 -c 6 -f
