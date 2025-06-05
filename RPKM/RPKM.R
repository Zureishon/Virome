# Load necessary libraries
library(data.table)  # For efficient data handling

# Read the count matrix from file
count_matrix <- fread("general/counts.txt", skip = 1, header = TRUE)  # Skip the metadata line

# Extract relevant columns (Geneid and counts for each sample)
gene_ids <- count_matrix$Geneid
counts <- as.matrix(count_matrix[, -(1:6)])  # Exclude non-count columns

####################################################Filter counts to keep only the ones with 10 or more reads mapped
###################This step does alter the data structure, but we are assuming thaat mapping less than 10 reads
#might be fishy, specially when others map more than 1000
counts[counts >= 1 & counts <= 10] <- 0

######Remove unnecessary part of the name

colnames(counts) <- sub("^co_sam/co_bam/", "", colnames(counts))

# Calculate total mapped reads per sample
total_reads <- colSums(counts)

# Calculate gene/transcript lengths (assuming Length column represents transcript length)
gene_lengths <- count_matrix$Length

# Calculate RPKM for each contig
rpkm_values <- matrix(NA, nrow = nrow(count_matrix), ncol = ncol(count_matrix) - 6)  # Initialize matrix for RPKM values

for (i in 1:nrow(count_matrix)) {
  for (j in 1:(ncol(count_matrix) - 6)) {
    reads <- counts[i, j]
    length <- gene_lengths[i]
    total <- total_reads[j]
    # Check for zero values to avoid NaN
    if (reads == 0 || total == 0) {
      rpkm_values[i, j] <- 0
    } else {
    # Calculate RPKM
    rpkm <- (reads / (length / 1000)) / (total / 1e6)
    
    # Assign RPKM value to matrix
    rpkm_values[i, j] <- rpkm
  }
  }
}

# Convert RPKM values to a data frame
rpkm_df <- as.data.frame(rpkm_values)
colnames(rpkm_df) <- colnames(counts)  # Assign sample names as column names
rpkm_df$Geneid <- gene_ids  # Add Geneid column

# Print or save RPKM data frame as needed
print(head(rpkm_df))  # Displaying the first few rows for inspection

# Example: Save RPKM values to a CSV file
write.csv(rpkm_df, file = "rpkm_values.csv", row.names = FALSE)


#########################################
############################################
###Now create the heatmap


# Load necessary libraries
library(pheatmap)

# Set row names to Gene IDs
rownames(rpkm_df) <- rpkm_df$Geneid
log2_rpkm_df <- rpkm_df[, -ncol(rpkm_df)]  # Remove Geneid column

# Transform the RPKM values to log2 scale
log2_rpkm_df <- log2(log2_rpkm_df + 1)  # Adding 1 to avoid log(0)

# Remove columns with constant values
log2_rpkm_df <- log2_rpkm_df[, apply(log2_rpkm_df, 2, function(x) sd(x, na.rm = TRUE) != 0)]

# Plot the heatmap
heatmap_result <- pheatmap(
  log2_rpkm_df,
  cluster_rows = TRUE,  # Clustering for rows
  cluster_cols = TRUE,  # Clustering for columns
  scale = "none",       # No scaling since already log-transformed
  clustering_distance_rows = "correlation",  # Use Pearson correlation for rows
  clustering_distance_cols = "correlation",  # Use Pearson correlation for columns
  clustering_method = "complete",
  color = colorRampPalette(c("white", "#F7F4F9", "#ff8969", "red"))(100),
  breaks = seq(0, 20, length.out = 101),
  main = "Heatmap of Log2 RPKM Values"
)


# Cut the dendrogram to form clusters
row_clusters <- cutree(heatmap_result$tree_row, k = 8)  # Example: create 11 clusters (adjust as desired)

# Prepare the annotation data frame
annotation_df <- data.frame(
  Cluster = factor(row_clusters, levels = 1:max(row_clusters))
)
rownames(annotation_df) <- rownames(log2_rpkm_df)

# Re-plot the heatmap with row annotations
pheatmap(
  log2_rpkm_df,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  clustering_distance_rows = "correlation",
  clustering_distance_cols = "correlation",
  clustering_method = "complete",
  color = colorRampPalette(c("white", "#F7F4F9", "#ff8969", "red"))(100),
  breaks = seq(0, 20, length.out = 101),
  annotation_row = annotation_df,
  main = "Heatmap of Log2 RPKM Values with Clustering"
)



###############################
################################ Get dataframe. Use those clusters for contig extension in the main script

# Create a dataframe with contigs and their assigned clusters
cluster_df <- data.frame(
  Contig = rownames(log2_rpkm_df),
  Cluster = factor(row_clusters, levels = 1:max(row_clusters))
)

# View the dataframe
print(cluster_df)

# Save the dataframe to a CSV file
write.csv(cluster_df, "contig_clusters.csv", row.names = FALSE)




#######################
## Now we make the contig clusters fasta files,
##Use those clusters of contigs for contig extension in the main script


library(Biostrings)

# Load the FASTA file
fasta_file <- readDNAStringSet("non_redundant_contigs.fasta")

# Create a list to store sequences by cluster
cluster_sequences <- split(names(fasta_file), cluster_df$Cluster)

# Iterate over each cluster and save the corresponding sequences
for (cluster in names(cluster_sequences)) {
  # Get the contigs for the current cluster
  contigs_in_cluster <- cluster_sequences[[cluster]]
  
  # Filter the sequences that match contigs in the current cluster
  matched_sequences <- fasta_file[grepl(paste(contigs_in_cluster, collapse = "|"), names(fasta_file))]
  
  # Save the sequences to a new FASTA file for each cluster
  writeXStringSet(matched_sequences, paste0("cluster_", cluster, ".fasta"))
}


###################################################


