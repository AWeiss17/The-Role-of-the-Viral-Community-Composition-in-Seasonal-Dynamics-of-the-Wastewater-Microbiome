########################################### Viral PERMANOVA ##########################################

####################################
### Load Required Libraries
#####################################

library(ggplot2)
library(vegan)
library(dplyr)
library(readr)
library(tibble)


#####################################
### Step 1: Load rarefied counttables & environmental tables
#####################################

data_rarefied_refseq <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_refseq_486.csv")

# Environmental Data
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")
Environmental_Data_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metatranscriptomics.csv")

# All PCoA results
PCoA_ARG <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_ARG.csv")

# PCoA result Eukaryotes Relativized!
PCoA_Euk_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Euk_rarefied_Rel.csv")

# PCoA result Prokaryotes k=2 relativized!!!!!
PCoA_Prok_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Prok_rarefied_Rel.csv")


# Set SRR_ID as column
data_rarefied_refseq <- data_rarefied_refseq %>%
  rename(SRR_ID = ...1)

Environmental_Data_DNA <- Environmental_Data_DNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten

# Remove outlier samples
removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
data_rarefied_refseq <- data_rarefied_refseq %>% filter(!SRR_ID %in% removed_refseq_samples)
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)

removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)


#removed_ARG_samples <- c("SRR1611149", "SRR1544596", "SRR1611146")
#Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_ARG_samples)
#data_rarefied_refseq <- data_rarefied_refseq %>% filter(!SRR_ID %in% removed_ARG_samples)

#removed_ARG_samples_RNA <- c("SRR1611150", "SRR1544599", "SRR1611147")
#Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_ARG_samples_RNA)
#PCoA_Prok_rarefied <- PCoA_Prok_rarefied %>% filter(!SRR_ID %in% removed_ARG_samples_RNA)
#PCoA_Euk_rarefied <- PCoA_Euk_rarefied %>% filter(!SRR_ID %in% removed_ARG_samples_RNA)


#####################################
### Step 2: Combine column with Environmental table 
#####################################
# rename PC1 coumn that it is clear it comes from Virus
PCoA_ARG <- PCoA_ARG %>%
  select(SRR_ID, PC1, PC2) %>%
  rename(ARG_PC1 = PC1, ARG_PC2 = PC2)

# combine PCoA result with Environmental table
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  left_join(PCoA_ARG, by = "SRR_ID")



# add Date to PCoA to combine PCoA and DNA rarefied data over Date and nor SRR_ID because this is different
PCoA_Euk_rarefied <- PCoA_Euk_rarefied %>%
  left_join(Environmental_Data_RNA %>% select(SRR_ID, Date), by = "SRR_ID")

# rename PC1 column that it is clear it comes from Virus
PCoA_Euk_rarefied <- PCoA_Euk_rarefied %>%
  select(Date, PC1, PC2) %>%
  rename(Euk_PC1 = PC1, Euk_PC2 = PC2)

# combine PCoA result with Environmental table
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  left_join(PCoA_Euk_rarefied, by = "Date")



# add Date to PCoA to combine PCoA and DNA rarefied data over Date and nor SRR_ID because this is different
PCoA_Prok_rarefied <- PCoA_Prok_rarefied %>%
  left_join(Environmental_Data_RNA %>% select(SRR_ID, Date), by = "SRR_ID")

# rename PC1 coumn that it is clear it comes from Virus
PCoA_Prok_rarefied <- PCoA_Prok_rarefied %>%
  select(Date, PC1, PC2) %>%
  rename(Prok_PC1 = PC1, Prok_PC2 = PC2)

# combine PCoA result with Environmental table
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  left_join(PCoA_Prok_rarefied, by = "Date")


#####################################
### Step 3: Relativize Data
#####################################
# only relativize the data if you DO NOT use bray in adonis2
# Data is already relativized because dataset comes from rarefraction!

# Relativierung der Sequenzdaten
#virus_data_matrix <- virus_data_matrix / rowSums(virus_data_matrix)  # Jede Zelle durch die Zeilensumme teilen

# Konvertiere zu Matrix für adonis2
virus_data_matrix <- data_rarefied_refseq %>%
  column_to_rownames(var = "SRR_ID") %>%  # Setze SRR_ID als Zeilennamen
  as.matrix()


#####################################
### Step 4: Prepare Environmental Table
#####################################
# Relevante Umweltvariablen auswählen
env_data <- Environmental_Data_DNA %>%
  filter(SRR_ID %in% rownames(virus_data_matrix)) # Nur passende SRR_IDs behalten
  
# Sortiere `env_data` exakt nach der Reihenfolge von `virus_data_matrix`
env_data <- env_data[match(rownames(virus_data_matrix), env_data$SRR_ID), ]

# Setze SRR_ID als Zeilennamen
env_data <- env_data %>%
  column_to_rownames(var = "SRR_ID")

# Überprüfe, ob die Zeilen exakt übereinstimmen
identical(rownames(virus_data_matrix), rownames(env_data))


#####################################
### Step 4: PERMANOVA with adonis2
#####################################
# by = terms tests each variable based on the order in the Model
# the order of the variables is important

adonis2(virus_data_matrix ~ ARG_PC1 * Prok_PC1 * Euk_PC1  * Season, 
        data = env_data, method = "bray", permutations = 999, by = "terms")







##########################################################################################################
#########################                  Prokaryote PERMANOVA                  #########################
##########################################################################################################


####################################
### Load Required Libraries
#####################################

library(ggplot2)
library(vegan)
library(dplyr)
library(readr)
library(tibble)

#####################################
### Step 1: Load rarefied counttables & environmental tables
#####################################

# Table RNA rarefied
data_rarefied_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_RNA.csv")

# Counttable (for Gernera identification)
count_data_pr2 <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_pr2.csv")
count_data_silva <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_silva.csv")

# Environmental RNA table
Environmental_Data_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metatranscriptomics.csv")
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")

# PCoA result table cmdscale
PCoA_Virus_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Virus_rarefied.csv")

# This is the new PCoA result without Mimiviridae and Iridoviridae
#PCoA_Virus_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Virus_new.csv")

# PCoA result Eukaryotes cmdscale
#PCoA_Euk_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Euk_rarefied.csv")

# PCoA result Eukaryotes Relativized!
PCoA_Euk_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Euk_rarefied_Rel.csv")

# All PCoA results
PCoA_ARG <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_ARG.csv")


#####################################
### Step 2: Extract Prokaryotes from rarefied table
#####################################

# Set SRR_ID as column
data_rarefied_RNA <- data_rarefied_RNA %>%
  rename(SRR_ID = ...1)

# Extract unique Prokaryotic genus names from count_data_silva
prokaryote_genera <- unique(count_data_silva$Genus)

# Keep only prokaryotic genera that actually exist in data_rarefied_RNA
prokaryote_genera_present <- prokaryote_genera[prokaryote_genera %in% colnames(data_rarefied_RNA)]

# Now select only prokaryotic taxa from `data_rarefied_RNA`
data_prokaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(prokaryote_genera_present))  


#####################################
### Step 3: Remove Outlier from Metadata
#####################################

# Remove outlier samples (are already exluded from rarefied data table, son only exclude them from Env. table)
removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)

# Remove outlier samples
removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)

# select only neccessary Metadata columns
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten

Environmental_Data_DNA <- Environmental_Data_DNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten



#####################################
### Step 4: Combine column with Environmental table 
#####################################
# add Date to PCoA to combine PCoA and RNA rarefied data over Date and nor SRR_ID because this is different
PCoA_Virus_rarefied <- PCoA_Virus_rarefied %>%
  left_join(Environmental_Data_DNA %>% select(SRR_ID, Date), by = "SRR_ID")

# rename PC1 coumn that it is clear it comes from Virus
PCoA_Virus_rarefied <- PCoA_Virus_rarefied %>%
  select(Date, PC1, PC2) %>%
  rename(Virus_PC1 = PC1, Virus_PC2 = PC2)

# combine PCoA result with Environmental table
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  left_join(PCoA_Virus_rarefied, by = "Date")



# rename PC1 coumn that it is clear it comes from Virus
PCoA_ARG <- PCoA_ARG %>%
  select(SRR_ID, PC1, PC2) %>%
  rename(ARG_PC1 = PC1, ARG_PC2 = PC2)

PCoA_ARG <- PCoA_ARG %>%
  left_join(Environmental_Data_DNA %>% select(SRR_ID, Date), by = "SRR_ID") %>%
  select(-SRR_ID)

# combine PCoA result with Environmental table
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  left_join(PCoA_ARG, by = "Date")



# add Date to PCoA to combine PCoA and DNA rarefied data over Date and nor SRR_ID because this is different
PCoA_Euk_rarefied <- PCoA_Euk_rarefied %>%
  left_join(Environmental_Data_RNA %>% select(SRR_ID, Date), by = "SRR_ID")

# rename PC1 coumn that it is clear it comes from Virus
PCoA_Euk_rarefied <- PCoA_Euk_rarefied %>%
  select(Date, PC1, PC2) %>%
  rename(Euk_PC1 = PC1, Euk_PC2 = PC2)

# combine PCoA result with Environmental table
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  left_join(PCoA_Euk_rarefied, by = "Date")




#####################################
### Step 5: Prepare PERMANOVA
#####################################
# prepare data matrix
# make matrix for adonis2
data_prokaryote_matrix <- data_prokaryote %>%
  column_to_rownames(var = "SRR_ID") %>%  
  as.matrix()

# relativize per row
data_prokaryote_rel <- data_prokaryote_matrix / rowSums(data_prokaryote_matrix)


# prepare Environmental table
# sort the order of SRR_IDs that they fit for both tables
env_data <- Environmental_Data_RNA[match(rownames(data_prokaryote_rel), Environmental_Data_RNA$SRR_ID), ]

# SRR_ID as rownames
env_data <- env_data %>%
  column_to_rownames(var = "SRR_ID")

# check rownames
identical(rownames(data_prokaryote_rel), rownames(env_data))


#####################################
### Step 6: PERMANOVA
#####################################
# adonis2 with Bray-Curtis (so no previous normalization of data needed)

## For our study we use this Permanova:
adonis2(data_prokaryote_rel ~ ARG_PC1 * Euk_PC1 * Virus_PC1 * Season ,
        data = env_data, method = "bray", permutations = 999, by = "terms")





##########################################################################################################
#########################                  Eukaryote PERMANOVA                  #########################
##########################################################################################################


####################################
### Load Required Libraries
#####################################

library(ggplot2)
library(vegan)
library(dplyr)
library(readr)
library(tibble)

#####################################
### Step 1: Load rarefied counttables & environmental tables
#####################################

# Table RNA rarefied
data_rarefied_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_RNA.csv")

# Counttable (for Gernera identification)
count_data_pr2 <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_pr2.csv")
count_data_silva <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_silva.csv")

# Environmental RNA table
Environmental_Data_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metatranscriptomics.csv")
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")

# PCoA result table cmdscale
PCoA_Virus_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Virus_rarefied.csv")

# PCoA result Prokaryotes k=2 relativized!
PCoA_Prok_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Prok_rarefied_Rel.csv")

# All PCoA results
PCoA_ARG <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_ARG.csv")


#####################################
### Step 2: Extract Eukaryotes from rarefied table
#####################################

# Set SRR_ID as column
data_rarefied_RNA <- data_rarefied_RNA %>%
  rename(SRR_ID = ...1)

# Extract unique Prokaryotic genus names from count_data_silva
eukaryote_genera <- unique(count_data_pr2$Genus)

# Keep only prokaryotic genera that actually exist in data_rarefied_RNA
eukaryote_genera_present <- eukaryote_genera[eukaryote_genera %in% colnames(data_rarefied_RNA)]

# Now select only prokaryotic taxa from `data_rarefied_RNA`
data_eukaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(eukaryote_genera_present))  


#####################################
### Step 3: Remove Outlier from Metadata
#####################################

# Remove outlier samples (are already exluded from rarefied data table, son only exclude them from Env. table)
removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)

# Remove outlier samples
removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)

# select only neccessary Metadata columns
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten

Environmental_Data_DNA <- Environmental_Data_DNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten



#####################################
### Step 4: Combine column with Environmental table 
#####################################
# add Date to PCoA to combine PCoA and RNA rarefied data over Date and nor SRR_ID because this is different
PCoA_Virus_rarefied <- PCoA_Virus_rarefied %>%
  left_join(Environmental_Data_DNA %>% select(SRR_ID, Date), by = "SRR_ID")

# rename PC1 coumn that it is clear it comes from Virus
PCoA_Virus_rarefied <- PCoA_Virus_rarefied %>%
  select(Date, PC1, PC2) %>%
  rename(Virus_PC1 = PC1, Virus_PC2 = PC2)

# combine PCoA result with Environmental table
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  left_join(PCoA_Virus_rarefied, by = "Date")



# rename PC1 coumn that it is clear it comes from Virus
PCoA_ARG <- PCoA_ARG %>%
  select(SRR_ID, PC1, PC2) %>%
  rename(ARG_PC1 = PC1, ARG_PC2 = PC2)

PCoA_ARG <- PCoA_ARG %>%
  left_join(Environmental_Data_DNA %>% select(SRR_ID, Date), by = "SRR_ID") %>%
  select(-SRR_ID)

# combine PCoA result with Environmental table
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  left_join(PCoA_ARG, by = "Date")



# add Date to PCoA to combine PCoA and DNA rarefied data over Date and nor SRR_ID because this is different
PCoA_Prok_rarefied <- PCoA_Prok_rarefied %>%
  left_join(Environmental_Data_RNA %>% select(SRR_ID, Date), by = "SRR_ID")

# rename PC1 coumn that it is clear it comes from Virus
PCoA_Prok_rarefied <- PCoA_Prok_rarefied %>%
  select(Date, PC1, PC2) %>%
  rename(Prok_PC1 = PC1, Prok_PC2 = PC2)

# combine PCoA result with Environmental table
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  left_join(PCoA_Prok_rarefied, by = "Date")




#####################################
### Step 5: Prepare PERMANOVA
#####################################
# prepare data matrix
# make matrix for adonis2
data_eukaryote_matrix <- data_eukaryote %>%
  column_to_rownames(var = "SRR_ID") %>%  
  as.matrix()

# relativize per row
data_eukaryote_rel <- data_eukaryote_matrix / rowSums(data_eukaryote_matrix)


# prepare Environmental table
# sort the order of SRR_IDs that they fit for both tables
env_data <- Environmental_Data_RNA[match(rownames(data_eukaryote_rel), Environmental_Data_RNA$SRR_ID), ]

# SRR_ID as rownames
env_data <- env_data %>%
  column_to_rownames(var = "SRR_ID")

# check rownames
identical(rownames(data_eukaryote_rel), rownames(env_data))


#####################################
### Step 6: PERMANOVA
#####################################
# adonis2 with Bray-Curtis (so no previous normalization of data needed)

## For our study we use this Permanova:
adonis2(data_eukaryote_rel ~ Virus_PC1 * ARG_PC1 * Prok_PC1 * Season ,
        data = env_data, method = "bray", permutations = 999, by = "terms")



##########################################################################################################
#########################                  ARG PERMANOVA                  #########################
##########################################################################################################
# ARG = Antibiotika Resistance Genes extracted form Metagenomics dataset

####################################
### Load Required Libraries
#####################################

library(ggplot2)
library(vegan)
library(dplyr)
library(readr)
library(tibble)

#####################################
### Step 1: Load ARG normalized table & PCoA result Virus/ Euk & Env. table
#####################################

# ARG from 16S normalized table
ARG_data <- read_csv("Master/Metatranscriptomics/Excel_lists/ARG/ARG_data_combined.csv")

# Environmental tables
Environmental_Data_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metatranscriptomics.csv")
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")

# PCoA result table cmdscale
PCoA_Virus_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Virus_rarefied.csv")

# This is the new PCoA result without Mimiviridae and Iridoviridae
#PCoA_Virus_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Virus_new.csv")

# PCoA result Eukaryotes Relativized!
PCoA_Euk_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Euk_rarefied_Rel.csv")

# PCoA result Prokaryotes k=2 relativized!
PCoA_Prok_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Prok_rarefied_Rel.csv")



#####################################
### Step 2: Prepare ARG table
#####################################
# long forma with SRR_ID as row
ARG_long <- ARG_data %>%
  pivot_longer(cols = -Family, names_to = "SRR_ID", values_to = "Count")

#  Wide format with Family as colum
ARG_transformed <- ARG_long %>%
  pivot_wider(names_from = Family, values_from = Count, values_fill = 0)


#####################################
### Step 3: Remove outliers
#####################################

# Remove outlier samples
removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)

# Remove outlier samples
removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)
ARG_transformed <- ARG_transformed %>% filter(!SRR_ID %in% removed_refseq_samples)

# select only neccessary Metadata columns
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten

Environmental_Data_DNA <- Environmental_Data_DNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten



#####################################
### Step 4: Combine column with Environmental table 
#####################################
# rename PC1 column that it is clear it comes from Virus
PCoA_Virus_rarefied <- PCoA_Virus_rarefied %>%
  select(SRR_ID, PC1, PC2) %>%
  rename(Virus_PC1 = PC1, Virus_PC2 = PC2)

# combine PCoA result with Environmental table
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  left_join(PCoA_Virus_rarefied, by = "SRR_ID")

# now for Eukaryotes
# add Date to PCoA to combine PCoA Euk and PCoA Virus over Date and not SRR_ID because this is different
PCoA_Euk_rarefied <- PCoA_Euk_rarefied %>%
  left_join(Environmental_Data_RNA %>% select(SRR_ID, Date), by = "SRR_ID")

# rename PC1 coumn that it is clear it comes from Virus
PCoA_Euk_rarefied <- PCoA_Euk_rarefied %>%
  select(Date, PC1, PC2) %>%
  rename(Euk_PC1 = PC1, Euk_PC2 = PC2)

# combine PCoA result with Environmental table
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  left_join(PCoA_Euk_rarefied, by = "Date")



#now for Prokaryotes
# add Date to PCoA to combine PCoA Euk and PCoA Virus over Date and not SRR_ID because this is different
PCoA_Prok_rarefied <- PCoA_Prok_rarefied %>%
  left_join(Environmental_Data_RNA %>% select(SRR_ID, Date), by = "SRR_ID")

# rename PC1 coumn that it is clear it comes from Virus
PCoA_Prok_rarefied <- PCoA_Prok_rarefied %>%
  select(Date, PC1, PC2) %>%
  rename(Prok_PC1 = PC1, Prok_PC2 = PC2)

# combine PCoA result with Environmental table
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  left_join(PCoA_Prok_rarefied, by = "Date")



#####################################
### Step 5: Prepare PERMANOVA
#####################################
# prepare data matrix
# Konvertiere zu Matrix für adonis2
# set SRR_ID as rowname
ARG_matrix <- ARG_transformed %>% 
  column_to_rownames("SRR_ID") %>%  
  as.matrix()


# prepare Environmental table
# sort the order of SRR_IDs that they fit for both tables
env_data <- Environmental_Data_DNA[match(rownames(ARG_matrix), Environmental_Data_DNA$SRR_ID), ]

# SRR_ID as rowname
env_data <- env_data %>%
  column_to_rownames(var = "SRR_ID")

# Check rownames
identical(rownames(ARG_matrix), rownames(env_data))


#####################################
### Step 6: PERMANOVA
#####################################
# Relativize data
ARG_matrix <- ARG_matrix / rowSums(ARG_matrix)

# adonis2 with Bray-Curtis (so no previous normalization of data needed)
adonis2(ARG_matrix ~  Prok_PC1 * Virus_PC1 * Euk_PC1 * Season, 
        data = env_data, method = "bray", permutations = 999, by = "terms")

