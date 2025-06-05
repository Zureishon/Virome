####IMPORTANT
##########################Use this code to get the taxonomy data from NCBI
library(rentrez)
library(dplyr)

# Function to fetch taxonomy info from NCBI
fetch_taxonomy_info <- function(ids) {
  results <- lapply(ids, function(id) {
    search <- entrez_search(db = "nucleotide", term = id)
    if (search$count > 0) {
      summary <- entrez_summary(db = "nucleotide", id = search$ids[1])
      taxid <- summary$taxid
      organism <- summary$organism     ##########New. Delete if doesnt work. Also make sure to delete organism from further things
      species <- summary$title
      return(data.frame(subject_id = id, taxid = taxid, species = species, organism = organism, stringsAsFactors = FALSE))
    } else {
      return(data.frame(subject_id = id, taxid = NA, species = NA, organism = NA, stringsAsFactors = FALSE))
    }
  })
  do.call(rbind, results)
}


# Main function to update BLAST results with taxonomy information
update_blast_results_with_taxonomy <- function(input_csv, output_csv) {
  # Read the data
  data <- read.csv(input_csv, header = FALSE)
  colnames(data) <- c("query_id", "subject_id", "identity", "alignment_length", 
                      "mismatches", "gap_opens", "q_start", "q_end", 
                      "s_start", "s_end", "evalue", "bit_score")
  
  # Get unique subject_ids
  unique_ids <- unique(data$subject_id)
  
  # Fetch taxonomy info
  taxonomy_info <- fetch_taxonomy_info(unique_ids)
  
  # Merge taxonomy info with original data
  if (is.null(taxonomy_info) || nrow(taxonomy_info) == 0) {
    # No taxonomy info found, assign NA to virus_category and retain subject_id
    data$virus_category <- NA_character_
  } else {
  
  # Merge taxonomy info with original data, also catalog them as virus or no virus
data <- data %>%
    left_join(taxonomy_info, by = "subject_id") %>%
    group_by(query_id) %>%
    slice(1) %>%
    mutate(virus_category = ifelse(grepl("virus|virion|viridae|phage", species, ignore.case = TRUE), "virus", "no virus"))
  }

  # Save the updated data
  write.csv(data, output_csv, row.names = FALSE)
}

##Run the function for the csv blast files, it adds taxonomy information, the category for virus or no_virus.

update_blast_results_with_taxonomy("KWNR.csv", "KWNR_S148_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P14.csv", "P14_S10_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P15.csv", "P15_S15_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P17.csv", "P17_S17_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P1.csv", "P1_S45_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P2.csv", "P2_S46_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P3.csv", "P3_S47_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P4.csv", "P4_S4_L002_R1_001_VBlast.csv")
update_blast_results_with_taxonomy("P5.csv", "P5_S11_L002_R1_001_VBlast.csv") ###
update_blast_results_with_taxonomy("P7.csv", "P7_S6_L002_R1_001_VBlast.csv") ###

update_blast_results_with_taxonomy("contigs_final.csv", "Contigs_finalVBlast.csv")


## Code to rename contigs. It works even if there is not an updated blast file (example P6, which has contigs, but no blast result)
rename_contigs <- function(fasta_file, blast_file = NULL, new_fasta, popu) {
  library(seqinr)
  library(dplyr)
  
  # Read the original FASTA file
  sequences <- read.fasta(fasta_file)
  
  # If blast_file is NULL, create an empty data frame for blast_results
  if (is.null(blast_file)) {
    blast_results <- data.frame(query_id = character(length(sequences)),
                                species = rep("NA", length(sequences)),
                                organism = rep("NA", length(sequences)),
                                virus_category = rep("NA", length(sequences)),
                                stringsAsFactors = FALSE)
  } else {
    # Read the Blast results CSV file
    blast_results <- read.csv(blast_file)
  }
  
  # Prepare the mapping between current and new names with Blast information
  mapping <- data.frame(
    current_name = names(sequences),
    species = rep(NA, length(sequences)),
    virus_category = rep(NA, length(sequences)),
    organism = rep(NA, length(sequences)),
    stringsAsFactors = FALSE
  )
  
  # Align blast_results with sequences using match
  match_indices <- match(mapping$current_name, blast_results$query_id)
  mapping$species <- blast_results$species[match_indices]
  mapping$virus_category <- blast_results$virus_category[match_indices]
  mapping$organism <- blast_results$organism[match_indices]
  
  # Generate new names including species and virus_category information
  for (i in 1:nrow(mapping)) {
    # Find the index of the current sequence in the sequences list
    sequence_index <- which(names(sequences) == mapping$current_name[i])
    
    # Extract contig size from the sequence
    seq_string <- paste(sequences[[sequence_index]], collapse = "")
    
    # Calculate the length of the concatenated sequence string
    contig_size <- nchar(seq_string)
    
    # Create new name in the format >P2_contig001_WNV_virus
    new_name <- paste0(popu, "_contig", sprintf("%03d", i), "_", contig_size, "_", 
                       mapping$species[i], "_", mapping$virus_category[i])
    
    # Assign new name to mapping dataframe
    mapping$new_name[i] <- new_name
  }
  
  # Select only current_name and new_name columns
  mapping <- mapping[, c("current_name", "new_name", "organism")]
  
  # Write the mapping file to disk with popu in the filename
  mapping_filename <- paste0(popu, "_mapping.txt")
  write.table(mapping, mapping_filename, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
  
  # Write the modified FASTA file with updated contig names
  fasta_output <- paste0(new_fasta, ".fasta")
  write.fasta(sequences, fasta_output, names = mapping$new_name)
}

##Run the function for each case. It will match the taxonomy entries with the contigs files, and change each contig name.
rename_contigs("KWNR_S148_L002_R1_001_200.fasta", "KWNR_S148_L002_R1_001_VBlast.csv", "KWNR_S148_L002_R1_001_renamed", "KWNR")
rename_contigs("P14_S10_L002_R1_001_200.fasta", "P14_S10_L002_R1_001_VBlast.csv", "P14_S10_L002_R1_001_renamed", "P14")
rename_contigs("P15_S15_L002_R1_001_200.fasta", "P15_S15_L002_R1_001_VBlast.csv", "P15_S15_L002_R1_001_renamed", "P15")
rename_contigs("P17_S17_L002_R1_001_200.fasta", "P17_S17_L002_R1_001_VBlast.csv", "P17_S17_L002_R1_001_renamed" , "P17")
rename_contigs("P1_S45_L002_R1_001_200.fasta", "P1_S45_L002_R1_001_VBlast.csv", "P1_S45_L002_R1_001_renamed" , "P1")
rename_contigs("P2_S46_L002_R1_001_200.fasta", "P2_S46_L002_R1_001_VBlast.csv", "P2_S46_L002_R1_001_renamed" , "P2")
rename_contigs("P3_S47_L002_R1_001_200.fasta", "P3_S47_L002_R1_001_VBlast.csv", "P3_S47_L002_R1_001_renamed" , "P3")
rename_contigs("P4_S4_L002_R1_001_200.fasta", "P4_S4_L002_R1_001_VBlast.csv", "P4_S4_L002_R1_001_renamed", "P4")
rename_contigs("P5_S11_L002_R1_001_200.fasta", "P5_S11_L002_R1_001_VBlast.csv", "P5_S11_L002_R1_001_renamed", "P5")
rename_contigs("P6_S5_L002_R1_001_200.fasta", new_fasta="P6_S5_L002_R1_001_renamed", popu="P6")
rename_contigs("P7_S6_L002_R1_001_200.fasta", "P7_S6_L002_R1_001_VBlast.csv", "P7_S6_L002_R1_001_renamed", "P7")

###Run for final contigs

rename_contigs("contigs_final.fasta", "Contigs_finalVBlast.csv", "contigs_final_renamed", "final")




