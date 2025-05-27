#####################################
### Load Required Libraries
#####################################
library(dplyr)
library(readr)

#####################################
### Step 1: Load Nils's Metatranscriptomic Data
#####################################

counttable_metatranscr01_genus_new <- read_delim("Master/Metatranscriptomics/Excel_lists/Metatranscriptomic Tables/counttable_metatranscr01_genus_new.csv", 
                                                 delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Standardize column names
colnames(counttable_metatranscr01_genus_new) <- c("SRR_ID", "db", "Supergroup", "Phylum", "Class", "Order", "Family", "Genus", "Count")

# Load table SRR_ID Metatranscriptomics with Date
SRA_Metatranscriptomics <- read_delim("Master/Metatranscriptomics/Excel_lists/Metatranscriptomic Tables/SRA_Metatranscriptomics.csv", 
                                      delim = ";", escape_double = FALSE, trim_ws = TRUE)


# Combine Nils's couttable with metadata
combined_Metatranscriptomic <- inner_join(counttable_metatranscr01_genus_new, SRA_Metatranscriptomics, by = "SRR_ID")

# save combined table
write.csv(combined_Metatranscriptomic, "combined_Metatranscriptomics.csv", row.names = FALSE)




#####################################
### Step 2: Load and Process Metagenomics Data
#####################################

combined_kraken2_report <- read_csv("Master/Metatranscriptomics/Excel_lists/Metagenomic Tables/combined_kraken2_report.csv")

# Load table with SRR_ID and Date
Sample_List_with_date <- read_delim("Master/Metatranscriptomics/Excel_lists/Metagenomic Tables/Sample_List_with_date.csv", 
                                    delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Add Date and Sample_ID to combined_kraken2_report
combined_kraken2_report_with_date <- combined_kraken2_report %>%
  left_join(Sample_List_with_date, by = "SRR_ID")

# Remove only the row where Rank == "U" and TaxID == 0
filtered_metagenomics <- combined_kraken2_report_with_date %>%
  filter(!(Rank == "U" & TaxID == 0))

# save combined table
write.csv(filtered_metagenomics, "combined_Metagenomics.csv", row.names = FALSE)


#####################################
### Step 3: Reorganize Metagenomics table (like Nils's table)
#####################################

# show unique Ranks in table
unique_ranks <- unique(filtered_metagenomics$Rank)
print(unique_ranks)
# assign ranks to Taxon
rank_map <- c(
  "R" = "Root",
  "D" = "Domain",
  "D1" = "Subdomain1",
  "D2" = "Subdomain2",
  "D3" = "Subdomain3",
  "D4" = "Subdomain4",
  "D5" = "Subdomain5",
  "K" = "Kingdom",
  "P" = "Phylum",
  "P1" = "Subphylum1",
  "C" = "Class",
  "C1" = "Subclass1",
  "O" = "Order",
  "O1"= "Suborder1",
  "F" = "Family",
  "F1" = "Subfamily1",
  "G" = "Genus",
  "G1" = "Subgenus1",
  "G2" = "Subgenus2",
  "S" = "Species",
  "S1" = "Subspecies1",
  "S2" = "Subspecies2",
  "S3" = "Subspecies3"
)

# add "Level"-Column based on Ranks
metagenomics <- filtered_metagenomics %>%
  mutate(Level = rank_map[Rank])

# create empty columns for each Level (after Rank-Mapping)
metagenomics_hierarchy <- metagenomics %>%
  mutate(
    Root = ifelse(Level == "Root", Taxon, NA),
    Domain = ifelse(Level == "Domain", Taxon, NA),
    Kingdom = ifelse(Level == "Kingdom", Taxon, NA),
    Phylum = ifelse(Level == "Phylum", Taxon, NA),
    Class = ifelse(Level == "Class", Taxon, NA),
    Family = ifelse(Level == "Family", Taxon, NA),
    Genus = ifelse(Level == "Genus", Taxon, NA),
    Species = ifelse(Level == "Species", Taxon, NA)
  ) %>%
  arrange(SRR_ID, Sample_ID, Date, match(Level, names(rank_map)))  # Sort for hierarchy

# Go through each row and fill after hierarchy
metagenomics_hierarchy <- metagenomics_hierarchy %>%
  group_by(SRR_ID, Sample_ID, Date) %>%
  arrange(SRR_ID, Sample_ID, Date) %>%  # sort rows within a group after hierarchy
  fill(Root, Domain, Kingdom, Phylum, Class, Family, Genus, Species, .direction = "down") %>%  # fill hierarchy columns to the bottom
  ungroup()

# add NA to Suffix (_X, _XX, _XXX, ...)
metagenomics_X <- metagenomics_hierarchy %>%
  mutate(
    Domain = ifelse(is.na(Domain), paste0("_X"), Domain),
    Kingdom = ifelse(is.na(Kingdom), paste0(Domain, "_X"), Kingdom),
    Phylum = ifelse(is.na(Phylum), paste0(Kingdom, "_X"), Phylum),
    Class = ifelse(is.na(Class), paste0(Phylum, "_X"), Class),
    Family = ifelse(is.na(Family), paste0(Class, "_X"), Family),
    Genus = ifelse(is.na(Genus), paste0(Family, "_X"), Genus),
    Species = ifelse(is.na(Species), paste0(Genus, "_X"), Species)
  )

# remove columns Percentage, Total_Reads, TaxID, Taxon, Level, Root
metagenomics_new <- metagenomics_X %>%
  select(-Percentage, -Total_Reads, -TaxID, -Taxon, -Level, -Root)

# add column "db" with "RefSeq"
metagenomics_db <- metagenomics_new %>%
  mutate(db = "RefSeq")

#delete first row
metagenomics_db <- metagenomics_db[-1, ]

# rename "Assigned_Reads" to "Count"
colnames(metagenomics_db)[which(colnames(metagenomics_db) == "Assigned_Reads")] <- "Count"


# Promote "uncultured" taxa to the next higher classification
cleaned_metagenomics <- metagenomics_db %>%
  mutate(
    Genus  = ifelse(Genus == "uncultured", paste0("uncultured ", Family), Genus),
    Family = ifelse(Family == "uncultured", paste0("uncultured ", Class), Family),
    Class  = ifelse(Class == "uncultured", paste0("uncultured ", Phylum), Class)
  )


# Create "Counttable" Metagenomics, Counts only for Family Level (filtered through Rank = F)
# Filter for Family-level assigned reads

family_assigned_reads <- cleaned_metagenomics %>%
  filter(Rank == "F")  # Family level only

# Remove column Rank, Species
metagenomics_new <- family_assigned_reads %>%
  select(-Rank, -Species)

metagenomics_final <- metagenomics_new
# save combined table
write.csv(metagenomics_final, "Metagenomics_final.csv", row.names = FALSE)




#####################################
### Step 4: Combine both Tables Metagenomics & Metatranscriptomics
#####################################

###### At first Column names and order need to be modified
########### Metatranscriptomics table ################
# Add column "Domain" and "Kingdom" before "Supergroup"
combined_Metatranscriptomic$Domain <- NA      
combined_Metatranscriptomic$Kingdom <- NA     

# Modify column order
Metatranscriptomic <- combined_Metatranscriptomic[, c("SRR_ID", "db", "Domain", "Kingdom", "Supergroup", 
                                                      "Phylum", "Class", "Order", "Family", "Genus", 
                                                      "Count", "Date", "Sample_ID")]

write.csv(Metatranscriptomic, "Metatranscriptomic.csv", row.names = FALSE)

########### Metagenomics table ################
# Add column "Supergroup" before "Phylum"
metagenomics_final$Supergroup <- NA      
# Add column "Order" before "Family"
metagenomics_final$Order <- NA

# Modify column order
Metagenomic <- metagenomics_final[, c("SRR_ID", "db", "Domain", "Kingdom", "Supergroup", 
                                      "Phylum", "Class", "Order", "Family", "Genus", 
                                      "Count", "Date", "Sample_ID")]

write.csv(Metagenomic, "Metagenomic.csv", row.names = FALSE)
