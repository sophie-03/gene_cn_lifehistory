# R script to run PGLS on cancer data and proteome orthogroup data

## FOR CLOUD CLUSTER

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Extract command-line arguments
col_name <- args[1]
output <- args[2]

setwd("/home/smatthews/placement/proteome_pgls")

library(phylolm)
library(readxl)
library(dplyr)
library(ape)
library(tidytree)
library(tibble)
library(tidyr)
library(geiger)
source("zach_pgls_function.R")
#source file: https://github.com/sophie-03/cancerAcrossVertebrates/blob/main/pglsPagel.R 

### ORTHOGROUP DATA
ortho <- read.table("Orthogroups_proteomes.tsv", fill = TRUE, sep = "\t", header = TRUE)
# Function to count genes in a cell
count_genes <- function(cell) {
  if (cell == "") {
    return(0)
  } else {
    genes <- unlist(strsplit(cell, ", "))
    return(length(genes))
  }
}
# Apply the function to each cell of the data frame
counts <- as.data.frame(lapply(ortho[, -1], function(col) sapply(col, count_genes)))
# Create the final output data frame
counts <- cbind(ortho$Orthogroup, counts)
#transpose df so species are rows and gene count is columns
counts <- t(counts)
# Set the column names
colnames(counts) <- counts[1, ]
# Remove the first row (which is now the column names)
counts <- counts[-1, ]
counts <- as.data.frame(counts)

#exclude species with a low busco score
#read in BUSCO results
busco <- read.table("busco_summaries.txt", header = TRUE)
#species to exclude
low_busco <- busco$Species[busco$C < 50]
counts <- counts[!(rownames(counts) %in% low_busco), ]

#this is where I removed low variance groups - but that was incorrect
filtered_data <- counts

### PHENOTYPE DATA
spreadsheet <- read.table("Mammals_Neoplasia_records.csv", header = TRUE, sep = ",")                            
#filter species to proteome species
filtered_spreadsheet <- spreadsheet[spreadsheet$Species %in% rownames(filtered_data), ]

#read in proteome sequence counts
proteome_size <- read.table("proteome_sequence_counts.txt", header = TRUE)
#remove .faa from species names
proteome_size$Species <- sub("\\..*", "", proteome_size$Species)
#attach proteome size to other data
filtered_spreadsheet <- left_join(filtered_spreadsheet, proteome_size, by = "Species")

#tree
tree <- read.tree("Mammals_tree_all_10.nwk")
#prune tree for just species we have data for
all_species <- tree$tip.label
data_species <- rownames(counts)
pruned_species <- setdiff(all_species, data_species)
pruned_tree <- drop.tip(tree, pruned_species)

### FORMAT DATA
subset_ortho <- filtered_data
# Convert row names to a new column called 'label'
subset_ortho <- data.frame(label = rownames(subset_ortho), subset_ortho, row.names = NULL)
#change Species column to label
filtered_spreadsheet <- filtered_spreadsheet %>%
  rename(label = Species)

## FILTER ORTHOGROUPS FOR ONES WITH A COPY NUMBER >0 FOR >50% SPECIES
# make all columns (except species column) numeric
#subset_ortho[, -1] <- lapply(subset_ortho[, -1], as.numeric)                               
# count how many species in each orthogroup are 0
#zero_counts <- colSums(subset_ortho[, -1] == 0)
# Find orthogroups where the number of zeros is less than 50% of the total species
#selected_orthogroups <- names(zero_counts[zero_counts < 0.5 * nrow(subset_ortho)])   
# Subset the dataframe for selected orthogroups
#subset_ortho <- subset_ortho[, c("label", selected_orthogroups)]                                


### PGLS: TRAIT AND ORTHOGROUPS
pgls_df <- data.frame(OG = character(), coeff = numeric(), r_squared = numeric(), lambda = numeric(), p_value = numeric(), stringsAsFactors = FALSE)


# Loop through each column of subset_ortho starting from the second column
set.seed(123)
for (col in names(subset_ortho)[-1]) {
tryCatch({

    # Initialize variables to keep track of outliers and removed species
    removed_species <- character(0)
    has_outliers <- TRUE

    # Subset data1 for the current column and age
    orthogroup <- subset_ortho[, c("label", col)]
    data1 <- left_join(orthogroup, filtered_spreadsheet, by = "label")
    data1[[col_name]] <- as.numeric(data1[[col_name]])
    data1[[col]] <- as.numeric(data1[[col]])
    rownames(data1) <- data1$label
    #SE
    SE<-setNames(data1$SE,data1$label)[rownames(data1)]
    
     #remove species with NA or negative values in the phenotype column
    data1 <- data1 %>%
      filter(!is.na(data1[[col_name]]) & data1[[col_name]] >= 0)
  
    while (has_outliers) {
      
      # pgls
      formula <- as.formula(paste0(col_name, " ~ ", col, " + sequence_count"))
      pgls <- pglsSEyPagel(formula, data=data1, tree=pruned_tree, se=SE, method = "ML")
      
      # Identify outliers in the pgls model
      res <- residuals(pgls, phylo = TRUE)
      res <- res / sqrt(var(res))[1]
      outliers <- abs(res) > 3
      
      if (any(outliers)) {
        # Add outliers to the list of removed species
        removed_species <- c(removed_species, names(outliers)[outliers])
        
        # Remove outliers from data1
        data1 <- data1[!rownames(data1) %in% removed_species, ]

        # pgls
        formula <- as.formula(paste0(col_name, " ~ ", col, " + sequence_count"))
        pgls <- pglsSEyPagel(formula, data=data1, tree=pruned_tree, se=SE, method = "ML")
      
        # Identify outliers in the pgls model
        res <- residuals(pgls, phylo = TRUE)
        res <- res / sqrt(var(res))[1]
        outliers <- abs(res) > 3
        
      } else {
        # No outliers found, exit the loop
        has_outliers <- FALSE
      }
    }
    
    # Extract and store statistics 
    p_value <- summary(pgls)$tTable[2,4]
    lambda <- summary(pgls)$modelStruct$corStruct
    std_error <- summary(pgls)$tTable[2,2]
    t_value <- summary(pgls)$tTable[2,3]
    coeff <- summary(pgls)$tTable[2,1]
    residual_stderror <- summary(pgls)$sigma

  #print statistics in output
    pgls_df <- rbind(pgls_df, data.frame(OG = col,
                                              coeff = coeff,
                                              p_value = p_value,
                                              lambda = lambda,
                                              std_error = std_error,
                                              t_value = t_value,
                                              residual_stderror = residual_stderror,
                                              outliers = paste(removed_species, collapse = ",")))
 }, error = function(e) {
    # Handle the error as needed, or simply continue to the next iteration
    cat("Error occurred in column:", col, "\n")
    print(e)
  })
}

#save output    note to self: the iteration of this script that has filters orthogroups, the output is written to with_min_species folder (so don't panic they're gonna overwrite)                          
output_filename <- paste("results/pgls_indv_", output, ".txt", sep = "")
write.table(pgls_df, output_filename, quote = FALSE, row.names = FALSE)
