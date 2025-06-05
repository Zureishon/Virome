# Load necessary libraries
library(data.table)  # For efficient data handling
library(dplyr)  # For clean data manipulation


# Read the count matrix from file
count_matrix <- fread("general/counts.txt", skip = 1, header = TRUE)  # Skip the metadata line

gene_id_mapping <- data.frame(
  original_id = c("MH188052.1", "NC_040716.1", "NC_032231.1", "MH188050.1", 
                  "MF176248.1", "MW434901.1", "MK628543.1"),
  new_id = c("Culex_Bunyavirus_2", "Culex_Iflavi-like_virus_4", "Hubei_mosquito_virus_4", 
             "Partitivirus-like_Culex_mosquito_virus", "Wuhan_Mosquito_Virus_6", 
             "Marma_virus", "Culex_narnavirus_1")
)

# Optional: Replace gene IDs using dplyr
count_matrix <- count_matrix %>%
  left_join(gene_id_mapping, by = c("Geneid" = "original_id")) %>%
  mutate(Geneid = ifelse(is.na(new_id), Geneid, new_id)) %>%  # Replace old Geneid with new_id
  select(-new_id)  # Remove the extra column

gene_ids <- count_matrix$Geneid

####################################Otherwise just use this code (if doesnt want to change gene ids)
###################################
# Extract relevant columns (Geneid and counts for each sample)
#gene_ids <- count_matrix$Geneid

counts <- as.matrix(count_matrix[, -(1:6)])  # Exclude non-count columns

####################################################Filter counts to keep only the ones with more than 300 reads mapped
###################This step does alter the data structure, but we are assuming that mapping less than 300 reads
#might be fishy, specially when others map more than 10000
counts[counts >= 1 & counts <= 300] <- 0

##Remove unnecessary names
colnames(counts) <- sub("^co_sam2/co_bam/(.*?)(_.*)?$", "\\1", colnames(counts))


###########change the column names for the codes.



# Calculate total mapped reads per sample
total_reads <- colSums(counts)

# Calculate gene/transcript lengths (assuming Length column represents transcript length)
gene_lengths <- count_matrix$Length

# Calculate RPKM for each gene/transcript
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
row_clusters <- cutree(heatmap_result$tree_row, k = 3)  # Example: create 11 clusters (adjust as desired)

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


################################################### Heatmap with STATES


collection_sites <- data.frame(
  Population = paste0("P", 1:17),  # Concatenate "P" with numbers 1 to 17
  Code = c("UT1", "UT2", "UT3", "WA1", "WA2", "WA3", "TX1", "TX2", "CO1", "CO2", "CO3", 
           "CA1", "CA2", "CA3", "CA4", "CA5", "CA6"),  # Manually assigned codes
  Location = c("Lake point", "Stansbury Park", "Stansbury Park", "Yakima/Grandview", 
               "Yakima/Grandview", "Yakima/Grandview", "Lubbock", "Lubbock", "Loveland", 
               "Loveland", "Fort Collins", "Cottonwood, Shasta County", "Anderson, Shasta County", 
               "Cottonwood, Shasta County", "Yolo County", "Sacramento County", "Yolo County"),
  State = c("UT", "UT", "UT", "WA", "WA", "WA", "TX", "TX", "CO", "CO", "CO", 
            "CA", "CA", "CA", "CA", "CA", "CA"),
  Latitude = c(40.669406, 40.651018, 40.647411, 46.193421, 46.224977, 46.205611, 
               33.462700, 33.590464, 40.403933, 40.401394, 40.529927, 40.414286, 
               40.445232, 40.396727, 38.610336, 38.107311, 38.810503),
  Longitude = c(-112.305263, -112.299081, -112.314970, -119.899400, -119.901053, 
                -119.892026, -101.891400, -101.938660, -105.071421, -105.108170, 
                -105.051447, -122.216945, -122.243092, -122.281994, -121.711942, 
                -121.649758, -121.808285)
)

# Print the updated data frame
print(collection_sites)

# Load necessary libraries
library(pheatmap)
library(dplyr)

# Map population to state
population_to_state <- collection_sites %>%
  select(Code, State) 

# Use colorblind-friendly colors for the states
state_colors <- c(
  "UT" = "#D55E00",  # Reddish-orange for Utah
  "WA" = "#0072B2",  # Blue for Washington
  "TX" = "#009E73",  # Bluish-green for Texas
  "CO" = "#F0E442",  # Yellow for Colorado
  "CA" = "#CC79A7"   # Pinkish purple for California
)

# Create a column annotation based on state
column_annotation <- data.frame(
  State = factor(population_to_state$State)
)
rownames(column_annotation) <- population_to_state$Code

# Define the annotation colors
annotation_colors <- list(State = state_colors)

##############Change names of the log2_rpkm_df

# Create a named vector for mapping
name_mapping <- setNames(collection_sites$Code, collection_sites$Population)

# Rename the columns of log2_rpkm_df
colnames(log2_rpkm_df) <- ifelse(colnames(log2_rpkm_df) %in% names(name_mapping), 
                                 name_mapping[colnames(log2_rpkm_df)], 
                                 colnames(log2_rpkm_df))

# Print the updated column names
colnames(log2_rpkm_df)


# Plot the heatmap with column annotation
heatmap_result <- pheatmap(
  log2_rpkm_df,
  cluster_rows = TRUE,   # Clustering for rows
  cluster_cols = TRUE,   # Clustering for columns
  scale = "none",        # No scaling since already log-transformed
  clustering_distance_rows = "correlation",   # Use Pearson correlation for rows
  clustering_distance_cols = "correlation",   # Use Pearson correlation for columns
  clustering_method = "complete",
  color = colorRampPalette(c("white", "#F7F4F9", "#ff8969", "red"))(100),
  breaks = seq(0, 20, length.out = 101),
  main = "Heatmap of Log2 RPKM Values",
  annotation_col = column_annotation,  # Add column annotation for states
  annotation_colors = annotation_colors
)


