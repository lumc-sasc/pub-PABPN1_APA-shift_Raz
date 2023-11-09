library(dplyr)      # Load the dplyr library for data manipulation
library(readxl)     # Load the readxl library for reading Excel files

# Read the Excel file into a data frame called INPUT
INPUT <- read_excel("C:/your/path/file.xlsx")

# Filter rows where 'feature' column contains "ENST"
# Then, filter columns that are not equal to 0 in all columns except the 'feature' column
# Set the 'feature' column as row names
INPUT <- INPUT %>%
  filter(grepl("ENST", feature)) %>%
  filter_at(vars(-feature), any_vars(. != 0)) %>%
  tibble::column_to_rownames("feature")

# Calculate the sums of reads in 3’UTR (proximal + distal) for each transcript, per genotype. 1:10 and 11:20 and create new columns CTRL and COI
INPUT$CTRL <- rowSums(INPUT[, 1:10])
INPUT$COI <- rowSums(INPUT[, 11:20])

# Filterout the rows where CTRL <10 and COI <10 (considering minimum 1 read/sample)
filtered_data <- INPUT[INPUT$CTRL > 10 & INPUT$COI > 10, ]

# The CTRL and COI rowSums columns are only used for filtering and can be removed
INPUT <- filtered_data[, -c(21, 22)]

# Select the proximal an distal columns 
prx <- INPUT[, c(1:5, 11:15)]
dist <- INPUT[, c(6:10, 16:20)]

# Create a binary vector 'y' to represent the groups (0 and 1)
y <- c(rep(0, 5), rep(1, 5))

# Calculate proximal to distal ratio in a log2 
APA_ratio <- as.matrix(log2((prx + 1) / (dist + 1)))

# Calculate p-values between CTRL and COI for each row (transcript) using a t-test
pvals <- apply(APA_ratio, 1, function(x) {t.test(x ~ y, paired = FALSE)$p.value})

# Adjust p-values for multiple comparisons using the False Discovery Rate (FDR) method
padj <- p.adjust(pvals, "fdr")

# Calculate log-fold change (logFC) for each row
logFC <- apply(APA_ratio, 1, function(x) diff((t.test(x ~ y)$est)))

# Combine logFC, p-values, and adjusted p-values into a data frame
d <- cbind(logFC, pvals, padj)
d <- as.data.frame(d)

# Write the results to an output CSV file
write.csv(d, file = "output.csv", row.names = TRUE)