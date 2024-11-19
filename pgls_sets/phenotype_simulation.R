# this script simulates a phenotype and runs pgls on this phenotype against the gene count of genes in set NR_TGFB

setwd("~/placement/geneset_pgls")

library(phylolm)
library(readxl)
library(dplyr)
library(ape)
library(tidytree)
library(tibble)
library(tidyr)
library(geiger)

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

col = "GOBP_NEGATIVE_REGULATION_OF_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION"
col_name="simulated"
                               
# Subset data1 for the current column and age
orthogroup <- geneset_counts[, c("label", col)]
data1 <- left_join(orthogroup, tibble_tree, by = "label")
data1[[col]] <- as.numeric(data1[[col]])
rownames(data1) <- data1$label
                               
# Number of simulations
n_simulations <- 1000
# Placeholder for p-values
p_values <- numeric(n_simulations)

# Loop through simulations
for (i in 1:n_simulations) {
  
  # Generate a random phenotype
  simulated_phenotype <- rTraitCont(pruned_tree)
  # Match simulated phenotype to data1
  data1$simulated <- simulated_phenotype[match(data1$label, names(simulated_phenotype))]
  
  # Perform PGLS with the simulated phenotype
  model <- phylolm(
    simulated ~ GOBP_NEGATIVE_REGULATION_OF_TRANSFORMING_GROWTH_FACTOR_BETA_PRODUCTION + mouse_equiv_genes,
    data = data1,
    phy = pruned_tree,
    model = "lambda"
  )
  
  # Extract the p-value for the first predictor
  p_values[i] <- summary(model)$coefficients[2,4]  # Adjust the index based on predictor order
}

write.table(p_values,"phenotype_simulation_pvalues.txt", row.names=FALSE, quote = FALSE, sep = "\t")
