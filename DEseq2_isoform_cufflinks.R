# Load necessary libraries
library(dplyr)
library(DESeq2)
library(stringr)

# Change the path to the actual path of your CSV file
counts_path <- "C:/Users/Qiang/Desktop/SIAT/data/mus_cufflinks_nobiascorrect_isoforms_combined_coverage.csv"

# Read in the data from the CSV file
count_data <- read.csv(counts_path, header = TRUE, sep = ",", check.names = FALSE)

# Inspect loaded data to confirm it is as expected
head(count_data)

print(nrow(count_data))
print(head(count_data))

# The input columns will be changed according to the different comparing groups
required_columns <- c("tracking_id", "ADULT1_coverage", "ADULT2_coverage", "ADULT3_coverage", "KO1_coverage", "KO2_coverage", "KO3_coverage")
missing_columns <- setdiff(required_columns, names(count_data))
if(length(missing_columns) > 0) {
  stop(paste("Missing required columns:", paste(missing_columns, collapse=", ")))
}

# Select the required columns for analysis
selected_count <- dplyr::select(count_data, tracking_id,
                                `ADULT1_coverage`,
                                `ADULT2_coverage`,
                                `ADULT3_coverage`,
                                `KO1_coverage`,
                                `KO2_coverage`,
                                `KO3_coverage`)

# Print the resulting data frame to check the selection
print(head(selected_count))

# Filtering the counts
# Exclude the first column (tracking_id column) when calculating row means
selected_count <- selected_count[rowMeans(selected_count[-1]) > 10,]

# Print the head of the selected_count data frame after filtering
print(nrow(selected_count))

# Ensure there are no NA values in the count matrix
count_matrix <- as.matrix(round(selected_count[,-1])) # Exclude the non-numeric 'tracking_id' column
rownames(count_matrix) <- selected_count$tracking_id  # Assuming 'tracking_id' are the gene/feature identifiers

# Remove rows with NA values
count_matrix <- na.omit(count_matrix)

# Create the sample information data frame directly in the script
sample_info <- data.frame(
  Sample = c("ADULT1_coverage", "ADULT2_coverage", "ADULT3_coverage", "KO1_coverage", "KO2_coverage", "KO3_coverage"),
  Type = c("musADULT", "musADULT", "musADULT", "musKO", "musKO", "musKO")
)
rownames(sample_info) <- sample_info$Sample
sample_info$Sample <- NULL

# Preview sample information to confirm proper load
print(head(sample_info))

# Set 'Type' as a factor
sample_info$Type <- as.factor(sample_info$Type)

# Check to ensure that the column names of 'count_matrix' match the row names of 'sample_info'
if (!all(rownames(sample_info) %in% colnames(count_matrix))) {
  mismatched_names <- setdiff(rownames(sample_info), colnames(count_matrix))
  print(mismatched_names)
  stop("Mismatch between sample names and selected count data columns")
}

# Create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = sample_info,
                              design = ~ Type)

# View dimensions of the DESeqDataSet
print(dim(dds))

# Filtering low count genes
dds <- dds[rowSums(counts(dds)) > 1,]
print(nrow(dds)) # Number of remaining features after filtering

# Run DESeq
dep <- DESeq(dds)

# Extracting results
res <- results(dep)
# Preview results
print(head(res))

# Filtering out results with NAs
diff <- na.omit(res)
print(dim(diff)) # Dimensions after filtering NAs

# Add 'tracking_id' as a new column to the results, match by the row names
diff$tracking_id <- rownames(diff)

# Write the final table to a CSV file, including 'tracking_id'
output_file <- "C:/Users/Qiang/Desktop/SIAT/data/DEseq2_musADULT_KO_cufflinks_isoform_result_counthreshold10.csv"
write.csv(diff, output_file, row.names = FALSE)
