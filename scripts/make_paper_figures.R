library(tidyverse)

############### functions ###################################################################################################
extract_transcript_id <- function(attributes) {
	match <- str_match(attributes, 'transcript_id "([^"]+)')
	if (!is.na(match[1, 2])) {
		return(match[1, 2])
	} else {
		return(NULL)
	}
}

# function to calculate relative transcript abundance
relative_transcript_abundance <- function(df) {
	df_sums <- df %>%
		select(gene_id, total_CPM) %>%
		group_by(gene_id) %>%
		summarize(total_CPM_gene = sum(total_CPM), .groups = "drop")

	merged_df <- df %>%
		inner_join(df_sums, by = "gene_id") %>%
		mutate(relative_abundance_percent = (total_CPM / total_CPM_gene) * 100,
			total_CPM_transcript = total_CPM) %>%
		select(-total_CPM)

  return(merged_df)
}

# function to fix column names
# TODO double check chatGPT code
fix_column_names <- function(df, is_gene = FALSE) {

	# Check if this is a gene counts object
	if (is_gene) {

		# Get count column names and create a list of new column names
		count_columns <- colnames(df)
		list_new_names <- c("gene_id")

		# Move gene_id from index to first column
		df <- df %>% rownames_to_column("gene_id") %>% select(-rowname)
	} else {
		# Set count columns and create a list of new names
		count_columns <- colnames(df)[3:length(colnames(df))]
		list_new_names <- c("transcript_id", "gene_id")
	}

	# Fix names one by one and add to the list of new names
	for (i in 1:length(count_columns)) {
		count_columns[i] <- gsub("_mapped", "_counts", count_columns[i])
		list_new_names <- c(list_new_names, count_columns[i])
	}

	# Rename columns
	colnames(df) <- list_new_names

	return(df)
}

# function to parse dataframe columns
parse_df_columns <- function(df, is_ref = TRUE, is_transcript = FALSE, is_prot = FALSE) {
  
	if (is_ref) {
		# Get gene ids
		df$gene_id <- str_match(df$other, "gene_id \"([^\"]+)\";")[2]
    
		# Get gene names
		df$gene_name <- str_match(df$other, "gene_name \"([^\"]+)\";")[2]
    
		# Get gene biotype
		df$gene_biotype <- str_match(df$other, "gene_biotype \"([^\"]+)\";")[2]
    
		if (is_transcript) {
			# Get transcript id and biotype
			df$transcript_id <- str_match(df$other, "transcript_id \"([^\"]+)\";")[2]
			df$transcript_biotype <- str_match(df$other, "transcript_biotype \"([^\"]+)\";")[2]
      
			if (is_prot) {
				# Get protein id, ccds id, and exon number
				df$protein_id <- str_match(df$other, "protein_id \"([^\"]+)\";")[2]
				df$ccds_id <- str_match(df$other, "ccds_id \"([^\"]+)\";")[2]
				df$exon_number <- str_match(df$other, "exon_number \"([^\"]+)\";")[2]
			}
		}
	} else {
		# Get gene id
		df$gene_id <- str_match(df$other, "gene_id \"([^\"]+)\";")[2]
    
		# Get transcript id
		df$transcript_id <- str_match(df$other, "transcript_id \"([^\"]+)\";")[2]
    
		# Get exon number
		df$exon_number <- str_match(df$other, "exon_number \"([^\"]+)\";")[2]
  
		# Drop "other", "dot_1", and "dot_2" columns
		df <- df %>% select(-c(other, dot_1, dot_2))
	}

	# Drop "other", "dot_1", and "dot_2" columns
	# Replace empty cells with NA
	df <- df 
		%>% select(-c(other, dot_1, dot_2)) 
		%>% mutate_all(~ifelse(. == "", NA, .))
  
	return(df)
}
############################################################################################################################
#################### Load Datasets ########################################################################################

column_names = c("seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attributes")

# 2014
# load the gtf, select only the transcript lines, and create a transcript_id column
df2014 <- read_tsv('../references/Homo_sapiens.GRCh38.76.gtf', col_names = column_names, comment = "#") %>%
	filter(feature == 'transcript') %>%
	mutate(transcript_id = extract_transcript_id(attributes)) %>%
	mutate(year = '2014')

# 2015
# load the gtf, select only the transcript lines, and create a transcript_id column
df2015 <- read_tsv('../references/Homo_sapiens.GRCh38.79.gtf', col_names = column_names, comment = "#") %>%
	filter(feature == 'transcript') %>%
	mutate(transcript_id = extract_transcript_id(attributes)) %>%
	mutate(year = '2015')

# 2016
# load the gtf, select only the transcript lines, and create a transcript_id column
df2016 <- read_tsv('../references/Homo_sapiens.GRCh38.85.gtf', col_names = column_names, comment = "#") %>%
	filter(feature == 'transcript') %>%
	mutate(transcript_id = extract_transcript_id(attributes)) %>%
	mutate(year = '2016')

# 2017
# load the gtf, select only the transcript lines, and create a transcript_id column
df2017 <- read_tsv('../references/Homo_sapiens.GRCh38.88.gtf', col_names = column_names, comment = "#") %>%
	filter(feature == 'transcript') %>%
	mutate(transcript_id = extract_transcript_id(attributes)) %>%
	mutate(year = '2017')

# 2018
# load the gtf, select only the transcript lines, and create a transcript_id column
df2018 <- read_tsv('../references/Homo_sapiens.GRCh38.92.gtf', col_names = column_names, comment = "#") %>%
	filter(feature == 'transcript') %>%
	mutate(transcript_id = extract_transcript_id(attributes)) %>%
	mutate(year = '2018')

# 2019
# load the gtf, select only the transcript lines, and create a transcript_id column
df2019 <- read_tsv('../references/Homo_sapiens.GRCh38.96.gtf', col_names = column_names, comment = "#") %>%
	filter(feature == 'transcript') %>%
	mutate(transcript_id = extract_transcript_id(attributes)) %>%
	mutate(year = '2019')

# 2020
# load the gtf, select only the transcript lines, and create a transcript_id column
df2020 <- read_tsv('../references/Homo_sapiens.GRCh38.99.gtf', col_names = column_names, comment = "#") %>%
	filter(feature == 'transcript') %>%
	mutate(transcript_id = extract_transcript_id(attributes)) %>%
	mutate(year = '2020')

# 2021
# load the gtf, select only the transcript lines, and create a transcript_id column
df2021 <- read_tsv('../references/Homo_sapiens.GRCh38.103.gtf', col_names = column_names, comment = "#") %>%
	filter(feature == 'transcript') %>%
	mutate(transcript_id = extract_transcript_id(attributes)) %>%
	mutate(year = '2021')

# 2022
# load the gtf, select only the transcript lines, and create a transcript_id column
df2022 <- read_tsv('../references/Homo_sapiens.GRCh38.106.gtf', col_names = column_names, comment = "#") %>%
	filter(feature == 'transcript') %>%
	mutate(transcript_id = extract_transcript_id(attributes)) %>%
	mutate(year = '2022')

# 2023
# load the gtf, select only the transcript lines, and create a transcript_id column
df2023 <- read_tsv('../references/Homo_sapiens.GRCh38.109.gtf', col_names = column_names, comment = "#") %>%
	filter(feature == 'transcript') %>%
	mutate(transcript_id = extract_transcript_id(attributes)) %>%
	mutate(year = '2023')

# ensure 2014 and 2015 only have the chr that match in 2018
df2014 <- df2014 %>%
	filter(seqname %in% df2018$seqname)
df2015 <- df2015 %>%
	filter(seqname %in% df2018$seqname)

##########################################################################################################################

# number of unique IDs
unique_transcript_counts <- tribble(
				   ~year, ~n_unique_transcripts,
				   '2014', length(unique(df2014 %>% pull(transcript_id))),
				   '2015', length(unique(df2015 %>% pull(transcript_id))),
				   '2016', length(unique(df2016 %>% pull(transcript_id))),
				   '2017', length(unique(df2017 %>% pull(transcript_id))),
				   '2018', length(unique(df2018 %>% pull(transcript_id))),
				   '2019', length(unique(df2019 %>% pull(transcript_id))),
				   '2020', length(unique(df2020 %>% pull(transcript_id))),
				   '2021', length(unique(df2021 %>% pull(transcript_id))),
				   '2022', length(unique(df2022 %>% pull(transcript_id))),
				   '2023', length(unique(df2023 %>% pull(transcript_id))) )
difference_between_years <- unique_transcript_counts %>%
	mutate(difference = n_unique_transcripts - lag(n_unique_transcripts)) %>%
	replace_na(difference = 0)

pdf("annotation_differences.pdf", oneFile = TRUE)

ggplot(unique_transcript_counts, aes(x = year, y = n_unique_transcripts, color = "#F8766D")) %>%
	geom_bar() +
	xlab('Year') +
	ylab('Number of Transcripts') +



