#############################################################
############################################################

library(GenomicRanges)
library(Rsamtools)
library(GenomicAlignments)
library(ggplot2)
library(reshape2)  # Load the reshape2 library for data reshaping

# Define the function. This will create the coverage files once we call it
plot_coverage_from_bam <- function(bam_file, desired_length) {
  
  # Strip the '.bam' extension to get the base name for the plot title and output file name
  base_name <- tools::file_path_sans_ext(bam_file)
  
  # Use BamFile object with index for efficient reading
  bam <- BamFile(bam_file)
  
  # Load alignments using ScanBamParam to work with the index
  param <- ScanBamParam(what=c("seq", "pos", "strand", "qwidth"))
  
  # Read alignments and filter by desired length
  bam_alignments <- readGAlignments(bam, param=param)
  bam_alignments <- bam_alignments[width(bam_alignments) %in% desired_length]
  
  # Separate reads by strand
  bam_positive <- bam_alignments[strand(bam_alignments) == "+"]
  bam_negative <- bam_alignments[strand(bam_alignments) == "-"]
  
  # Calculate coverage for each strand
  coverage_positive <- coverage(bam_positive)
  coverage_negative <- coverage(bam_negative)
  
  # Check the number of contigs with non-zero coverage
  non_zero_pos_cov <- sum(sapply(coverage_positive, function(x) sum(x > 0)))
  non_zero_neg_cov <- sum(sapply(coverage_negative, function(x) sum(x > 0)))
  print(paste("Number of contigs with non-zero positive coverage:", non_zero_pos_cov))
  print(paste("Number of contigs with non-zero negative coverage:", non_zero_neg_cov))
  
  # Extract coverage values for contigs with non-zero reads
  pos_cov_values <- unlist(coverage_positive[sapply(coverage_positive, function(x) sum(x > 0)) > 0])
  neg_cov_values <- unlist(coverage_negative[sapply(coverage_negative, function(x) sum(x > 0)) > 0])
  
  # Check if there are any positive or negative coverage values
  print(paste("Number of non-zero positive coverage values:", length(pos_cov_values)))
  print(paste("Number of non-zero negative coverage values:", length(neg_cov_values)))
  
  # Pad the shorter vector with zeros to match the length and convert to numeric
  max_length <- max(length(pos_cov_values), length(neg_cov_values))
  pos_cov_values <- as.numeric(c(pos_cov_values, rep(0, max_length - length(pos_cov_values))))
  neg_cov_values <- as.numeric(c(neg_cov_values, rep(0, max_length - length(neg_cov_values))))
  
  # Create a data frame for plotting
  df <- data.frame(
    Position = seq_along(pos_cov_values),
    Positive = pos_cov_values,
    Negative = -neg_cov_values  # This is now valid
  )
  
  # Reshape the data frame for ggplot2
  df_long <- melt(df, id.vars = "Position")
  
  # Generate the plot
  p <- ggplot(df_long, aes(x = Position, y = value, fill = variable)) +
    geom_area(alpha = 1, position = 'identity') +  # Smooth filled area
    scale_fill_manual(values = c("#6baed6", "#fb6a4a")) +  # Use the preferred light blue/red colors
    geom_hline(yintercept = 0, color = "grey") +  # Solid line at coverage = 0
    labs(title = paste(base_name), x = "Position", y = "Coverage") +
    theme_minimal() +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      panel.border = element_blank(),      # Remove panel border
      axis.line = element_line(),          # Keep axis lines (brackets)
      axis.ticks = element_line(color = "black"),  # Add axis ticks (indents)
      axis.ticks.length = unit(0.15, "cm")  # Length of the ticks (indentations)
    )
  
  # Save the plot as PNG
  output_file <- paste0(base_name, "_length_", desired_length, ".pdf") ####################Change to png if needed
  ggsave(output_file, plot = p, width = 4, height = 2, device = cairo_pdf)   ################change as needed, dpi = 100
  
  print(paste0("Plot saved as ", output_file))
}

# Calling the function requires the presence of the .bam and .bai (it makes it faster) files. Example usage:

#Calls for each virus with read length 21

plot_coverage_from_bam("Culex_Bunyavirus_2.bam", 21)
plot_coverage_from_bam("Culex_Iflavi-like_virus_4.bam", 21)
plot_coverage_from_bam("Hubei_mosquito_virus_4.bam", 21)
plot_coverage_from_bam("Partitivirus-like_Culex_mosquito_virus.bam", 21)
plot_coverage_from_bam("Wuhan_Mosquito_Virus_6.bam", 21)
plot_coverage_from_bam("Marma_virus.bam", 21)
plot_coverage_from_bam("Culex_narnavirus_1.bam", 21)

#Calls for each virus with read lengths 24 to 29
plot_coverage_from_bam("Culex_Bunyavirus_2.bam", 24:29)
plot_coverage_from_bam("Culex_Iflavi-like_virus_4.bam", 24:29)
plot_coverage_from_bam("Hubei_mosquito_virus_4.bam", 24:29)
plot_coverage_from_bam("Partitivirus-like_Culex_mosquito_virus.bam", 24:29)
plot_coverage_from_bam("Wuhan_Mosquito_Virus_6.bam", 24:29)
plot_coverage_from_bam("Marma_virus.bam", 24:29)
plot_coverage_from_bam("Culex_narnavirus_1.bam", 24:29)


###############################################
##################################################

