
# this script uses the number of sequences in the proteome that ALSO have an equivalent mouse gene (and are therefore included in the genesets) as a covariate. 
# this is instead of using the raw sequence count as a covariate which is done in pgls_sets_reproductive_traits.R

setwd("~/placement/geneset_pgls")

library(phylolm)
library(readxl)
library(dplyr)
library(ape)
library(tidytree)
library(tibble)
library(tidyr)
library(geiger)

# Retrieve command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Extract command-line arguments
col_name <- args[1]
output <- args[2]

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


#data
spreadsheet <- read_excel("Mammals_Final_List.xlsx")
                               
#read in proteome sequence counts
proteome_size <- read.table("proteome_sequence_counts.txt", header = TRUE)
#remove .faa from species names
proteome_size$Species <- sub("\\..*", "", proteome_size$Species)
#attach proteome size to other data
spreadsheet <- left_join(spreadsheet, proteome_size, by = "Species")

## make birthweight/adultweight column
# Replace negative values with NA in birth_weight.g. and adult_weight.g.
spreadsheet <- spreadsheet %>%
  mutate(birth_weight.g. = ifelse(birth_weight.g. < 0, NA, birth_weight.g.),
         adult_weight.g. = ifelse(adult_weight.g. < 0, NA, adult_weight.g.))                            
#create column for the birth weight to adult weight ratio
spreadsheet$birth_adult_weight_ratio <- spreadsheet$birth_weight.g. / spreadsheet$adult_weight.g.

##make log(adult weight) column
spreadsheet$log_adult_weight <- log(spreadsheet$adult_weight.g.)
                               
#tree
tree <- read.tree("Mammals_tree_all_10.nwk")
#prune tree for just species we have data for
all_species <- tree$tip.label
data_species <- rownames(counts)
pruned_species <- setdiff(all_species, data_species)
pruned_tree <- drop.tip(tree, pruned_species)
#make a tree tibble
tibble_tree <- as_tibble(pruned_tree)

## add data to the tree
#read tree
data <- spreadsheet[spreadsheet$Species %in% pruned_tree$tip.label, ]
#rename Species column to label
colnames(data)[colnames(data) == "Species"] <- "label"
#add data to tree tibble
tibble_tree <- left_join(tibble_tree, data, by = "label")
#convert to treedata object
treedat <- as.treedata(tibble_tree)

#rename Species to label
colnames(geneset_counts)[colnames(geneset_counts) == "Species"] <- "label"


# Create an empty dataframe to store p-values
pgls_df <- data.frame(OG = character(), coeff = numeric(), r_squared = numeric(), lambda = numeric(), p_value = numeric(), stringsAsFactors = FALSE)

# Loop through each column of subset_ortho starting from the second column
set.seed(123)
for (col in names(geneset_counts)[-1]) {
  tryCatch({
    cat("Processing column:", col, "\n")
    # Initialize variables to keep track of outliers and removed species
    removed_species <- character(0)
    has_outliers <- TRUE
    
    # Subset data1 for the current column and age
    orthogroup <- geneset_counts[, c("label", col)]
    data1 <- left_join(orthogroup, tibble_tree, by = "label")
    data1[[col_name]] <- as.numeric(data1[[col_name]])
    data1[[col]] <- as.numeric(data1[[col]])
    rownames(data1) <- data1$label
    
    #remove species with NA or negative values in the phenotype column
    data1 <- data1 %>%
      filter(!is.na(data1[[col_name]]) & data1[[col_name]] >= 0)
    
    while (has_outliers) {
      
      # Run phylolm for the current column and age
      pgls <- phylolm(data1[[col_name]] ~ data1[[col]] + mouse_equiv_genes + max_longevity.months., data1, phy = pruned_tree, model = "lambda")
      
      # Identify outliers in the pgls model
      res <- residuals(pgls, phylo = TRUE)
      res <- res / sqrt(var(res))[1]
      outliers <- abs(res) > 3
      
      if (any(outliers)) {
        # Add outliers to the list of removed species
        removed_species <- c(removed_species, names(outliers)[outliers])
        
        # Remove outliers from data1
        data1 <- data1[!rownames(data1) %in% removed_species, ]

        # Run phylolm for the current column and age
      pgls <- phylolm(data1[[col_name]] ~ data1[[col]] + mouse_equiv_genes + max_longevity.months., data1, phy = pruned_tree, model = "lambda")
      
      # Identify outliers in the pgls model
      res <- residuals(pgls, phylo = TRUE)
      res <- res / sqrt(var(res))[1]
      outliers <- abs(res) > 3
        
      } else {
        # No outliers found, exit the loop
        has_outliers <- FALSE
      }
    }
    
    # Extract and store p-value in the p_values dataframe for lambda model
    p_value <- summary(pgls)$coefficients[2,4]
    r_squared <- pgls$r.squared
    lambda <- pgls$optpar
    coeff <- summary(pgls)$coefficients[2,1]
    std_error <- summary(pgls)$coefficients[2,2]
    
    # Store results in pgls_df dataframe
    pgls_df <- rbind(pgls_df, data.frame(OG = col,
                                          coeff = coeff,
                                          r_squared = r_squared,
                                          lambda = lambda,
                                          p_value = p_value,
                                          std_error = std_error,
                                          outliers = paste(removed_species, collapse = ",")))
  }, error = function(e) {
    # Handle the error as needed, or simply continue to the next iteration
    cat("Error occurred in column:", col, "\n")
    print(e)
  })
}


output_filename <- paste("results/covar_is_mouse_equiv_sequence_count/pgls_genesets_", output, ".txt", sep = "")

write.table(pgls_df, output_filename, quote = FALSE, row.names = FALSE)
