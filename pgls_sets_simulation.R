
# this script uses the number of sequences in the proteome that ALSO have an equivalent mouse gene (and are therefore included in the genesets) as a covariate. 


## FOR CLOUD CLUSTER

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Extract command-line arguments
col_name <- args[1]
output <- args[2]

setwd("/home/smatthews/placement/geneset_pgls")

library(phylolm)
library(readxl)
library(dplyr)
library(ape)
library(tidytree)
library(tibble)
library(tidyr)
library(geiger)
source("zach_pgls_function.R")

#read data
ortho <- read.table("OrthogroupsNov09.tsv", fill = TRUE, sep = "\t", header = TRUE)

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

mh_counts <- read.table("mh_set_counts.txt", header = TRUE)
m2cp_counts <- read.table("m2cp_set_counts.txt", header = TRUE)
m5_counts <- read.table("m5_set_counts.txt", header = TRUE)
geneset_counts <- merge(merge(mh_counts, m2cp_counts, by = "Species"), m5_counts, by = "Species") 


### PHENOTYPE DATA
spreadsheet <- read.table("Mammals_Neoplasia_records.csv", header = TRUE, sep = ",")                            


#read in proteome sequence counts
proteome_size <- read.table("proteome_sequence_counts.txt", header = TRUE)
#remove .faa from species names
proteome_size$Species <- sub("\\..*", "", proteome_size$Species)
#attach proteome size to other data
spreadsheet <- left_join(spreadsheet, proteome_size, by = "Species")

#tree
tree <- read.tree("Mammals_tree_all_10.nwk")
#prune tree for just species we have data for
all_species <- tree$tip.label
data_species <- rownames(counts)
pruned_species <- setdiff(all_species, data_species)
pruned_tree <- drop.tip(tree, pruned_species)
tibble_tree <- as_tibble(pruned_tree)

## add data to the tree
#read tree
data <- spreadsheet[spreadsheet$Species %in% pruned_tree$tip.label, ]
#rename Species column to label
colnames(data)[colnames(data) == "Species"] <- "label"
#add data to tree tibble
tibble_tree <- left_join(tibble_tree, data, by = "label")

#rename Species to label
colnames(geneset_counts)[colnames(geneset_counts) == "Species"] <- "label"

### FORMAT DATA
subset_ortho <- counts
# Convert row names to a new column called 'label'
subset_ortho <- data.frame(label = rownames(subset_ortho), subset_ortho, row.names = NULL)

## FILTER ORTHOGROUPS FOR ONES WITH A COPY NUMBER >0 FOR >50% SPECIES
# make all columns (except species column) numeric
subset_ortho[, -1] <- lapply(subset_ortho[, -1], as.numeric)                               
# count how many species in each orthogroup are 0
zero_counts <- colSums(subset_ortho[, -1] == 0)
# Find orthogroups where the number of zeros is less than 50% of the total species
selected_orthogroups <- names(zero_counts[zero_counts < 0.5 * nrow(subset_ortho)])   
# Subset the dataframe for selected orthogroups
subset_ortho <- subset_ortho[, c("label", selected_orthogroups)]   
                               
#VARIANCE
#read in variance of each orthogroup
variance <- read.table("cn_variability.txt")  
colnames(variance) <- c("Orthogroup", "variance")                               
# merge  orthogroups with gene names
mouse_ortho <- read.table("mouse_gene.orthogroup_data.tsv", header = TRUE, fill = TRUE)  
variance <- left_join(variance, mouse_ortho, by = "Orthogroup")                               
# remove orthogroups from variance that don't have >50% species with a non-zero copy number                               
variance <- variance %>% filter(Orthogroup %in% colnames(subset_ortho))

## READ IN GENE SETS AND THEIR GENES
mouse_m5 <- read.table("m5.all.v2023.2.Mm.symbols.gmt", header = FALSE, fill = TRUE, col.names = 1:10000)
#remove empty columns
mouse_m5 <- mouse_m5[, colSums(!is.na(mouse_m5)) > 0]
#remove column with web link
mouse_m5 <- mouse_m5[, -2]
#transpose
mouse_m5 <- t(mouse_m5)
mouse_m5 <- as.data.frame(mouse_m5)
#make set name the column name
colnames(mouse_m5) <- mouse_m5[1,]
#remove name row from df
mouse_m5 <- mouse_m5[-1,]
                               
## get variances of genes in set
#get genes in set
gene_list <- mouse_m5$GOBP_NEGATIVE_REGULATION_OF_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION
#subset variance list for genes of interest
gene_variance <- subset(variance, gene %in% gene_list)
#get min and max variance for this set
min_var <- min(gene_variance$variance)
max_var <- max(gene_variance$variance)         

 #subset genes for all genes with a variance between the min and the max
subset_variance <- variance[variance$variance >= min_var & variance$variance <= max_var, ]
#remove the genes that are in the gene set from this list
subset_variance <- subset_variance[!(subset_variance$gene %in% gene_list), ]
#remove orthogroup repeats
subset_variance <- subset_variance %>% distinct(Orthogroup, .keep_all = TRUE)                               
                               
### PGLS: TRAIT AND ORTHOGROUPS
pgls_df <- data.frame(OG = character(), coeff = numeric(), r_squared = numeric(), lambda = numeric(), p_value = numeric(), stringsAsFactors = FALSE)


# Loop through each column of subset_ortho starting from the second column
for (i in 1:1000) {
tryCatch({

    # Initialize variables to keep track of outliers and removed species
    removed_species <- character(0)
    has_outliers <- TRUE

     # Randomly sample 12 genes from the subset
  sampled_genes <- sample(unique(subset_variance$gene), 12)
  # Get the orthogroups of these genes
  sampled_variance <- subset_variance[subset_variance$gene %in% sampled_genes, ]
  # Subset count data for just these orthogroups
  subset_counts <- counts[, unique(sampled_variance$Orthogroup), drop = FALSE]
  # Numeric
  subset_counts <- as.data.frame(lapply(subset_counts, as.numeric))
  # Sum the counts for each species
  subset_counts$total <- rowSums(subset_counts)
  # Add total to a new dataframe
  totals <- as.data.frame(rownames(counts))
  colnames(totals) <- "label"
  totals$total <- subset_counts$total
  # Join total to phenotype data
  data2 <- left_join(totals, tibble_tree, by = "label")
  # Set species to row names
  rownames(data2) <- data2$label
  # SE
  SE <- setNames(data2$SE, data2$label)[rownames(data2)]
    
     #remove species with NA or negative values in the phenotype column
    data2 <- data2 %>%
      filter(!is.na(data2[[col_name]]) & data2[[col_name]] >= 0)
  
    while (has_outliers) {
      
      # pgls
      formula <- as.formula(paste0(col_name, " ~ total + mouse_equiv_genes"))
      pgls <- pglsSEyPagel(formula, data=data2, tree=pruned_tree, se=SE, method = "ML")
      
      # Identify outliers in the pgls model
      res <- residuals(pgls, phylo = TRUE)
      res <- res / sqrt(var(res))[1]
      outliers <- abs(res) > 3
      
      if (any(outliers)) {
        # Add outliers to the list of removed species
        removed_species <- c(removed_species, names(outliers)[outliers])
        
        # Remove outliers from data1
        data2 <- data2[!rownames(data2) %in% removed_species, ]

        # pgls
        formula <- as.formula(paste0(col_name, " ~ total + mouse_equiv_genes"))
        pgls <- pglsSEyPagel(formula, data=data2, tree=pruned_tree, se=SE, method = "ML")
      
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
    pgls_df <- rbind(pgls_df, data.frame(OG = i,
                                              coeff = coeff,
                                              p_value = p_value,
                                              lambda = lambda,
                                              std_error = std_error,
                                              t_value = t_value,
                                              residual_stderror = residual_stderror,
                                              outliers = paste(removed_species, collapse = ",")))
 }, error = function(e) {
    # Handle the error as needed, or simply continue to the next iteration
    cat("Error occurred in column:", i, "\n")
    print(e)
  })
}

#save output                              
output_filename <- paste("results/", output, ".txt", sep = "")
write.table(pgls_df, output_filename, quote = FALSE, row.names = FALSE)
