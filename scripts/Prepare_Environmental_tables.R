#####################################
### Load Required Libraries
#####################################
library(dplyr)
library(readr)


#####################################
### Step 1: Load Environmental data table
#####################################
Environmental_Data_raw <- read_delim("Master/Metatranscriptomics/Excel_lists/Environmental_Data_raw.csv", 
                                     delim = ";", escape_double = FALSE, trim_ws = TRUE)


# Table contains addition Metadata information about Nitrate and Phosphate
Metadata_NO3_PO4_Zusatz <- read_delim("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Metadata_NO3_PO4_Zusatz.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)


#####################################
### Step 2: Prepare Environmental data table
#####################################

Environmental_Data_raw$Date <- as.Date(Environmental_Data_raw$Date, format = "%d.%m.%Y")
Metadata_NO3_PO4_Zusatz$Date <- as.Date(Metadata_NO3_PO4_Zusatz$Date, format = "%d.%m.%Y")

# use Left JOIN for combinig both tables over common Date
Environmental_Data_raw <- Environmental_Data_raw %>%
  left_join(Metadata_NO3_PO4_Zusatz, by = "Date")

# add column with Month
Environmental_Data_raw$Month <- as.numeric(format(Environmental_Data_raw$Date, "%m"))

# Add a season column based on the month
Environmental_Data_raw$Season <- NA  # Create an empty column

# Assign seasons based on month ranges
Environmental_Data_raw$Season[Environmental_Data_raw$Month %in% c(3, 4, 5)] <- "spring"
Environmental_Data_raw$Season[Environmental_Data_raw$Month %in% c(6, 7, 8)] <- "summer"
Environmental_Data_raw$Season[Environmental_Data_raw$Month %in% c(9, 10, 11)] <- "autumn"
Environmental_Data_raw$Season[Environmental_Data_raw$Month %in% c(12, 1, 2)] <- "winter"

# Convert the Season column to a factor with ordered levels
Environmental_Data_raw$Season <- factor(Environmental_Data_raw$Season, levels = c("spring", "summer", "autumn", "winter"))

write.csv(Environmental_Data_raw, "Environmental_Data.csv", row.names = FALSE)

# load SRR_ID Metatrascriptomic to match with Environmental_Data (over SRR_ID)
SRA_Metatranscriptomics <- read_delim("Master/Metatranscriptomics/Excel_lists/Metatranscriptomic Tables/SRA_Metatranscriptomics.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)

# combine SRR_ID Metatrascriptomic with Environmental_Data
SRA_Metatranscriptomics$Date <- as.Date(SRA_Metatranscriptomics$Date, format = "%d.%m.%Y")
Environmental_Data_RNA <- merge(Environmental_Data_raw, SRA_Metatranscriptomics, by = "Date", all = TRUE)

write.csv(Environmental_Data_RNA, "Environmental_Data_Metatranscriptomics.csv", row.names = FALSE)


# load SRR_ID Metagenomic to match with Environmental_Data (over SRR_ID)
Sample_List_with_date <- read_delim("Master/Metatranscriptomics/Excel_lists/Metagenomic Tables/Sample_List_with_date.csv", 
                                    delim = ";", escape_double = FALSE, trim_ws = TRUE)

# combine SRR_ID Metagenomic with Environmental_Data
Sample_List_with_date$Date <- as.Date(Sample_List_with_date$Date, format = "%d.%m.%Y")
Environmental_Data_DNA <- merge(Environmental_Data_raw, Sample_List_with_date, by = "Date", all = TRUE)

write.csv(Environmental_Data_DNA, "Environmental_Data_Metagenomics.csv", row.names = FALSE)
