#######################################################
# Libraries
library(dplyr)
library(readr)
library(stringr)
library(data.table)
setwd("C:\\Users\\mekar\\Desktop\\Extracting structural features from viral proteins\\data") 


####################################################### Blast section
# Blast processing data
# Load BLAST CSV
blast <- read_csv("blast_results.csv")

# First set of filtering
blast_clean <- blast %>%
  mutate(
    Query_Cover_num = as.numeric(str_remove(`Query Cover`, "%")),
    Per_ident_num = as.numeric(`Per. ident`),
    E_value_num = as.numeric(`E value`),
    Acc_Len_num = as.numeric(`Acc. Len`),
    description_lower = str_to_lower(Description),
    is_partial = str_detect(description_lower, "partial"),
    is_relevant_protein = str_detect(description_lower, "coat|capsid|cp|3b protein")
  ) %>%
  filter( Query_Cover_num >= 85, E_value_num < 1e-5, Acc_Len_num >= 150, Acc_Len_num <= 230,
    is_partial == FALSE, is_relevant_protein == TRUE, !(is_partial == TRUE & Query_Cover_num < 80)
    ) %>%
  mutate( identity_category = case_when( Per_ident_num >= 90 ~ "very_close", Per_ident_num >= 70 ~ "close",
      Per_ident_num >= 50 ~ "medium", TRUE ~ "distant")
  )

# Keep best hit per virus
best_per_virus <- blast_clean %>%
  group_by(`Scientific Name`) %>%
  arrange( E_value_num, desc(Query_Cover_num), desc(Per_ident_num), desc(`Total Score`)) %>% slice(1) %>%
  ungroup()

# Selecting up to 10 representatives across identity categories
selected <- best_per_virus %>%
  group_by(identity_category) %>%
  arrange(E_value_num, desc(Query_Cover_num), desc(Per_ident_num)) %>%
  slice_head(n = 3) %>%
  ungroup() %>%
  arrange(desc(Per_ident_num)) %>%
  slice_head(n = 10)

selected_5 <- selected %>% select(Description,`Scientific Name`,E_value_num,Query_Cover_num,Per_ident_num,`Total Score`,Acc_Len_num, identity_category)
colnames(selected_5) <- c("description","scientific_name","e_value","q_cover","per_ident","total_score","acc_len","identity_category")

write_csv(selected_5, "C:\\Users\\mekar\\Desktop\\Extracting structural features from viral proteins\\outputs\\selected_5_blast_representatives.csv")

#######################################################
####################################################### Clustal omega section
# Read Clustal alignment file counting and separating the amino acids into exact, similar weakly comparisons
#in order to check were those amino acids are in the structure of the reference protein 
lines <- readLines("clustalo-I20260512-224057-0895-15794156-p2m.aln-clustal_num")

#Lines preprocessing --> remobing the CLUSTAL title line& empty lines
lines <- lines[!str_detect(lines, "^CLUSTAL")]
lines <- lines[lines != ""]

# Sequence IDs included in the alignment
seq_ids <- c( "YP_233104.1", "NP_689395.1", "UNZ99878.1","1CWP_A", "QXV86510.1")

# Storing the complete aligned sequence for each ID
aligned_sequences <- setNames(rep("", length(seq_ids)), seq_ids)
sequence_pattern <- paste0("^(", paste(seq_ids, collapse = "|"), ")")

# Reconstructing the full aligned sequences from the alignment blocks
for (line in lines) {
  if (str_detect(line, sequence_pattern)) {
    parts <- str_split(str_trim(line), "\\s+")[[1]]
    id <- parts[1]
    seq_part <- parts[2]
    aligned_sequences[id] <- paste0(aligned_sequences[id],seq_part)
  }
}

# Identifying where the aligned sequence begins & measuring the length of each alignment block
first_seq_line <- lines[ str_detect(lines, paste0("^", seq_ids[1])) ][1]
first_seq_part <- str_split( str_trim(first_seq_line),"\\s+")[[1]][2]
start_pos <- str_locate(first_seq_line,fixed(first_seq_part))[1, "start"]
sequence_block_lengths <- c()
for (line in lines[
  str_detect(lines, paste0("^", seq_ids[1]))]) {
  
  parts <- str_split(str_trim(line), "\\s+")[[1]]
  sequence_block_lengths <- c( sequence_block_lengths,nchar(parts[2]))
}

# Extracting the consensus symbols under each alignment block
consensus_clean <- c()
block_index <- 1

for (line in lines) {
  is_sequence_line <- str_detect( line, sequence_pattern)
  is_consensus_line <- ( !is_sequence_line && !str_detect(line, "^CLUSTAL") && line != "")
  if (
    is_consensus_line &&
    block_index <= length(sequence_block_lengths)
  ) {
    block_len <- sequence_block_lengths[block_index]

    if (nchar(line) < start_pos + block_len - 1) {
      
      line <- str_pad( line, width = start_pos + block_len - 1, side = "right")
    }
    
    consensus_piece <- substr( line, start_pos, start_pos + block_len - 1)
    consensus_clean <- c( consensus_clean, consensus_piece)
    block_index <- block_index + 1
  }
}

#Combining the extracted consensus blocks into one complete consensus sequence
consensus <- paste0( consensus_clean, collapse = "")

# Reference sequence, split into individual characters
ref_seq <- aligned_sequences["1CWP_A"] 
ref_chars <- strsplit( ref_seq,"")[[1]]
consensus_chars <- strsplit(consensus,"")[[1]]

# Mapping the alignment positions to the reference residue numbering while excluding alignment gaps
ref_position <- cumsum(ref_chars != "-")
ref_position[ ref_chars == "-"] <- NA

# Building a sequence-conservation table
conservation_table <- data.frame(alignment_position = seq_along(consensus_chars),
  reference_position_1CWP = ref_position, reference_residue = ref_chars,
  conservation_symbol = consensus_chars)

# Keeping only alignment positions that map to the residues in 1CWP
conservation_table <- conservation_table %>%
  filter(!is.na(reference_position_1CWP))

# Exact conserved positions
exact_conserved <- conservation_table %>%
  filter( conservation_symbol == "*")
strongly_similar <- conservation_table %>%
  filter( conservation_symbol == ":")
weakly_similar <- conservation_table %>%
  filter( conservation_symbol == ".")
conserved_or_similar <- conservation_table %>%
  filter( conservation_symbol %in% c("*", ":", "."))

# Save outputs
write.csv(conservation_table, "C:\\Users\\mekar\\Desktop\\Extracting structural features from viral proteins\\outputs\\all_alignment_positions_mapped_to_1CWP.csv", row.names = FALSE)
write.csv(exact_conserved, "C:\\Users\\mekar\\Desktop\\Extracting structural features from viral proteins\\outputs\\exact_conserved_positions_star.csv", row.names = FALSE)
write.csv(strongly_similar, "C:\\Users\\mekar\\Desktop\\Extracting structural features from viral proteins\\outputs\\strongly_similar_positions_colon.csv", row.names = FALSE)
write.csv(conserved_or_similar, "conserved_or_similar_positions.csv", row.names = FALSE)
write.csv(weakly_similar, "C:\\Users\\mekar\\Desktop\\Extracting structural features from viral proteins\\outputs\\weakly_similar_positions.csv", row.names = FALSE)


#########################################
######################################### Conservation Chimera "Strings"
# Making the residue strings for chain A to use them into ChimeraX for coloring each conserved residue on top of the chain
# Visual purposes
make_residue_string <- function(file) {
  df <- read_csv(file)
  residues <- df %>% filter(!is.na(reference_position_1CWP)) %>%
    pull(reference_position_1CWP) %>% unique() %>% sort()
  paste(residues, collapse = ",")
}

exact_residues <- make_residue_string("C:\\Users\\mekar\\Desktop\\Extracting structural features from viral proteins\\outputs\\exact_conserved_positions_star.csv")
strong_residues <- make_residue_string("C:\\Users\\mekar\\Desktop\\Extracting structural features from viral proteins\\outputs\\strongly_similar_positions_colon.csv")
weak_residues <- make_residue_string("C:\\Users\\mekar\\Desktop\\Extracting structural features from viral proteins\\outputs\\weakly_similar_positions.csv")


#Copying them for Chimera
cat(paste0("/A:", exact_residues, "\n\n"))
cat(paste0("/A:", strong_residues, "\n\n"))
cat(paste0("/A:", weak_residues, "\n\n"))


######################################### Dali section
#########################################
read_dali_file <- function(file, query_chain, database) {
  df <- read_csv(file, show_col_types = FALSE)
  df %>% mutate( query_chain = query_chain, database = database)
  }


dali_C <- read_table("1cwpC.dali.tsv", comment = "#", col_names = TRUE, 
                     show_col_types = FALSE)

#Checking column names
print(colnames(dali_C))
head(dali_C)

# Cleaning & standardizing columns & Removing exact self-hit
dali_C_clean <- dali_C %>%
  rename(
    query = query,
    database = database,
    subject = sbjct,
    Z_score = `z-score`,
    rmsd = rmsd,
    ali_length = `ali-length`,
    subject_length = `sbjct-length`,
    seq_identity = `seq-identity`,
    subject_sequence = `sbjct-sequence`
  ) %>%
  mutate(
    subject_pdb = str_sub(subject, 1, 4),
    subject_chain = str_sub(subject, 5, 5)
  ) %>%
  filter(!subject %in% c("1cwpB", "1cwpC","1cwpA"))


# First Filtering Dali
dali_C_filtered <- dali_C_clean %>%
  filter( Z_score >= 8, rmsd <= 4, ali_length >= 80)

# Removing redundancy -->  keeping best hit per PDB entry (avoiding 
# keeping many chains from the same PDB structure)
dali_C_best_per_pdb <- dali_C_filtered %>%
  group_by(subject_pdb) %>%
  arrange(desc(Z_score), rmsd, desc(ali_length), desc(seq_identity)) %>%
  slice(1) %>%
  ungroup()

# Creating similarity categories like in blast 
dali_C_best_per_pdb <- dali_C_best_per_pdb %>%
  mutate(
    structural_similarity_category = case_when(
      Z_score >= 20 ~ "very_strong",
      Z_score >= 8  ~ "strong",
      TRUE ~ "weak"
    ),
    sequence_identity_category = case_when(
      seq_identity >= 90 ~ "very_close_sequence",
      seq_identity >= 70 ~ "close_sequence",
      seq_identity >= 50 ~ "medium_sequence",
      TRUE ~ "low_sequence"
    )
  )

#Selecting top representative hits (10)
selected_dali_C <- dali_C_best_per_pdb %>%
  arrange(desc(Z_score), rmsd, desc(ali_length)) %>%
  slice_head(n = 11) %>% select(everything(),-`compressed-sbjct-dssp`,-description,-description,-qstarts,-sstarts,
                                -lengths,-rotation,-translation,-all_pfam,-subject_pdb,-database,-query )

selected_dali_C <- selected_dali_C %>% select(everything(),-structural_similarity_category,-sequence_identity_category)

write_csv(selected_dali_C, "C:\\Users\\mekar\\Desktop\\Extracting structural features from viral proteins\\outputs\\structure_representatives.csv")
