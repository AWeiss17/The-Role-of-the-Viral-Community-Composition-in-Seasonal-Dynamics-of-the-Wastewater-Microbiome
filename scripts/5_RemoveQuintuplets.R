#####################################
### Load Required Libraries
#####################################
library(tidyr)
library(dplyr)
library(tibble)
library(readr)


#####################################
### Step 1: save counttables Metatranscriptomics from Nils created after renamer.R as csv
#####################################

txt_data <- read.table("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_pr2.txt", header = TRUE, sep = "", stringsAsFactors = FALSE, fill = TRUE)
write.csv(txt_data, "counttable_metatranscr01_genus_new_renamed_pr2.csv", row.names = FALSE)

txt_data <- read.table("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_silva.txt", header = TRUE, sep = "", stringsAsFactors = FALSE, fill = TRUE)
write.csv(txt_data, "counttable_metatranscr01_genus_new_renamed_silva.csv", row.names = FALSE)




#####################################
### Step 2: read counttables Metatranscriptomics & Metagenomics and do some renamings
#####################################

data_pr2 <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_pr2.csv")

data_silva <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_silva.csv")

data_refseq <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Metagenomic.csv")


# rename "accession" to "SRR_ID" in Nils counttables (pr2 & silva)
colnames(data_pr2)[which(colnames(data_pr2) == "accession")] <- "SRR_ID"
colnames(data_silva)[which(colnames(data_silva) == "accession")] <- "SRR_ID"

# rename "count" to "Count" in Nils counttables (pr2 & silva)
colnames(data_pr2)[which(colnames(data_pr2) == "count")] <- "Count"
colnames(data_silva)[which(colnames(data_silva) == "count")] <- "Count"

# Add Missing Columns to Nils's Count Tables
data_pr2$Domain <- NA
data_pr2$Kingdom <- NA

data_silva$Domain <- NA
data_silva$Kingdom <- NA


# Load table SRR_ID Metatranscriptomics with Date
SRA_Metatranscriptomics <- read_delim("Master/Metatranscriptomics/Excel_lists/Metatranscriptomic Tables/SRA_Metatranscriptomics.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Combine Metadata with Nils's Count Table
count_data_pr2 <- inner_join(data_pr2, SRA_Metatranscriptomics, by = "SRR_ID")
count_data_silva <- inner_join(data_silva, SRA_Metatranscriptomics, by = "SRR_ID")

# Reorder Columns (Ensure Consistency for Combination)
count_data_pr2 <- count_data_pr2[, c(
  "SRR_ID", "db", "Domain", "Kingdom", "Supergroup", "Phylum", "Class", 
  "Order", "Family", "Genus", "Count", "n", "Date", "Sample_ID"
)]

count_data_silva <- count_data_silva[, c(
  "SRR_ID", "db", "Domain", "Kingdom", "Supergroup", "Phylum", "Class", 
  "Order", "Family", "Genus", "Count", "n", "Date", "Sample_ID"
)]

count_data_pr2$Date <- as.Date(count_data_pr2$Date, format = "%d.%m.%Y")
count_data_silva$Date <- as.Date(count_data_silva$Date, format = "%d.%m.%Y")

# Save Modified Count Tables
write.csv(count_data_pr2, "count_data_pr2.csv", row.names = FALSE)
write.csv(count_data_silva, "count_data_silva.csv", row.names = FALSE)



#####################################
### Step 3: if there are homologous family names they are combined and the counts merged
#####################################
count_data_refseq <- data_refseq %>%
  group_by(SRR_ID, db, Domain, Kingdom, Supergroup, Phylum, Class, Order, Family, Genus, Date, Sample_ID) %>%
  dplyr::summarize(
    count = sum(Count),  # Sum up the count column for the grouping
    n = dplyr::n(),      # Add the n column
    .groups = "drop"
  )

colnames(count_data_refseq)[which(colnames(count_data_refseq) == "count")] <- "Count"
count_data_refseq$Date <- as.Date(count_data_refseq$Date, format = "%d.%m.%Y")

write.csv(count_data_refseq, "count_data_refseq.csv", row.names = FALSE)



#####################################
### Step 4: change to wide format
#####################################

## Rownames as SRR_ID
# PR2: spread function spreads out values (Count)
count_data_pr2_wide_SRR <- count_data_pr2 %>%
  select(SRR_ID, db, Genus, Count) %>%  # Select relevant columns
  pivot_wider(
    names_from = Genus,                 # Make Genera the column names
    values_from = Count,                # Fill columns with Counts
    values_fill = 0                     # Replace missing values with 0
  )

# Convert to data frame and set row names
count_data_pr2_wide_SRR <- as.data.frame(count_data_pr2_wide_SRR)  # Convert tibble to data frame
rownames(count_data_pr2_wide_SRR) <- count_data_pr2_wide_SRR$SRR_ID  # Set row names
count_data_pr2_wide_SRR <- count_data_pr2_wide_SRR[,-c(1,2)]  # Remove the SRR_ID column

# Filter out low counts in a column (e.g., singletons to quintupletons)
boxplot(log10(colSums(count_data_pr2_wide_SRR)), main = "Column Sums (Log10)")  # Visualize distribution
lowhitcount <- 5  # Threshold for exclusion

count_data_pr2_wide_SRR <- count_data_pr2_wide_SRR[, colSums(count_data_pr2_wide_SRR, na.rm = TRUE) > lowhitcount]

write.csv(count_data_pr2_wide_SRR, "count_data_pr2_wide_SRR.csv", row.names = TRUE)

# SILVA: spread function spreads out values (Count)
count_data_silva_wide_SRR <- count_data_silva %>%
  select(SRR_ID, db, Genus, Count) %>%  # Select relevant columns
  pivot_wider(
    names_from = Genus,                 # Make Genera the column names
    values_from = Count,                # Fill columns with Counts
    values_fill = 0                     # Replace missing values with 0
  )

# Convert to data frame and set row names
count_data_silva_wide_SRR <- as.data.frame(count_data_silva_wide_SRR)  # Convert tibble to data frame
rownames(count_data_silva_wide_SRR) <- count_data_silva_wide_SRR$SRR_ID  # Set row names
count_data_silva_wide_SRR <- count_data_silva_wide_SRR[,-c(1,2)]  # Remove the SRR_ID column

# Filter out low counts in a column (e.g., singletons to quintupletons)
boxplot(log10(colSums(count_data_silva_wide_SRR)), main = "Column Sums (Log10)")  # Visualize distribution
lowhitcount <- 5  # Threshold for exclusion

count_data_silva_wide_SRR <- count_data_silva_wide_SRR[, colSums(count_data_silva_wide_SRR, na.rm = TRUE) > lowhitcount]

write.csv(count_data_silva_wide_SRR, "count_data_silva_wide_SRR.csv", row.names = TRUE)

# RefSeq: spread function spreads out values (Count)
count_data_refseq_wide_SRR <- count_data_refseq %>%
  select(SRR_ID, db, Family, Count) %>%  # Select relevant columns
  pivot_wider(
    names_from = Family,                 # Make Genera the column names
    values_from = Count,                # Fill columns with Counts
    values_fill = 0                     # Replace missing values with 0
  )

# Convert to data frame and set row names
count_data_refseq_wide_SRR <- as.data.frame(count_data_refseq_wide_SRR)  # Convert tibble to data frame
rownames(count_data_refseq_wide_SRR) <- count_data_refseq_wide_SRR$SRR_ID  # Set row names
count_data_refseq_wide_SRR <- count_data_refseq_wide_SRR[,-c(1,2)]  # Remove the SRR_ID column

# Filter out low counts in a column (e.g., singletons to quintupletons)
boxplot(log10(colSums(count_data_refseq_wide_SRR)), main = "Column Sums (Log10)")  # Visualize distribution
lowhitcount <- 5  # Threshold for exclusion

count_data_refseq_wide_SRR <- count_data_refseq_wide_SRR[, colSums(count_data_refseq_wide_SRR, na.rm = TRUE) > lowhitcount]

write.csv(count_data_refseq_wide_SRR, "count_data_refseq_wide_SRR.csv", row.names = TRUE)



## Rownames as Dates
# PR2: spread function spreads out values (Count)
count_data_pr2_wide_Date <- count_data_pr2 %>%
  select(Date, db, Genus, Count) %>%  # Select relevant columns
  pivot_wider(
    names_from = Genus,                 # Make Genera the column names
    values_from = Count,                # Fill columns with Counts
    values_fill = 0                     # Replace missing values with 0
  )

# Convert to data frame and set row names
count_data_pr2_wide_Date <- as.data.frame(count_data_pr2_wide_Date)  # Convert tibble to data frame
rownames(count_data_pr2_wide_Date) <- count_data_pr2_wide_Date$Date  # Set row names
count_data_pr2_wide_Date <- count_data_pr2_wide_Date[,-c(1,2)]  # Remove the SRR_ID column

write.csv(count_data_pr2_wide_Date, "count_data_pr2_wide_Date.csv", row.names = TRUE)

# SILVA: spread function spreads out values (Count)
count_data_silva_wide_Date <- count_data_silva %>%
  select(Date, db, Genus, Count) %>%  # Select relevant columns
  pivot_wider(
    names_from = Genus,                 # Make Genera the column names
    values_from = Count,                # Fill columns with Counts
    values_fill = 0                     # Replace missing values with 0
  )

# Convert to data frame and set row names
count_data_silva_wide_Date <- as.data.frame(count_data_silva_wide_Date)  # Convert tibble to data frame
rownames(count_data_silva_wide_Date) <- count_data_silva_wide_Date$Date  # Set row names
count_data_silva_wide_Date <- count_data_silva_wide_Date[,-c(1,2)]  # Remove the SRR_ID column

write.csv(count_data_silva_wide_Date, "count_data_silva_wide_Date.csv", row.names = TRUE)

# RefSeq: spread function spreads out values (Count)
count_data_refseq_wide_Date <- count_data_refseq %>%
  select(Date, db, Family, Count) %>%  # Select relevant columns
  pivot_wider(
    names_from = Family,                 # Make Genera the column names
    values_from = Count,                # Fill columns with Counts
    values_fill = 0                     # Replace missing values with 0
  )

# Convert to data frame and set row names
count_data_refseq_wide_Date <- as.data.frame(count_data_refseq_wide_Date)  # Convert tibble to data frame
rownames(count_data_refseq_wide_Date) <- count_data_refseq_wide_Date$Date  # Set row names
count_data_refseq_wide_Date <- count_data_refseq_wide_Date[,-c(1,2)]  # Remove the SRR_ID column

write.csv(count_data_refseq_wide_Date, "count_data_refseq_wide_Date.csv", row.names = TRUE)




#####################################
### Step 5: combine all 3 tables
#####################################

# Ensure row names are explicitly stored as a column
count_data_pr2_wide_SRR <- tibble::rownames_to_column(count_data_pr2_wide_SRR, var = "SRR_ID")
count_data_silva_wide_SRR <- tibble::rownames_to_column(count_data_silva_wide_SRR, var = "SRR_ID")
count_data_refseq_wide_SRR <- tibble::rownames_to_column(count_data_refseq_wide_SRR, var = "SRR_ID")

# Combine PR2 and SILVA tables
count_data_combined <- full_join(count_data_pr2_wide_SRR, count_data_silva_wide_SRR, by = "SRR_ID")

# Combine with RefSeq table
count_data_wide <- full_join(count_data_combined, count_data_refseq_wide_SRR, by = "SRR_ID")

# Convert back to row names and remove extra column
rownames(count_data_wide) <- count_data_wide$SRR_ID
count_data_wide$SRR_ID <- NULL


# Save the final combined table
write.csv(count_data_wide, paste0("counttable_wide.csv"), row.names = TRUE, quote = FALSE)
