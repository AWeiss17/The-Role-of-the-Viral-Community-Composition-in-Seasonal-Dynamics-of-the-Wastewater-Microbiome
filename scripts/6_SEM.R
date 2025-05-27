#########################################################################################################################################
###########################################    Structural Equation Model (SEM) with Shannon  ##########################################
##################################################################################################################################

####################################
### Load Required Libraries
#####################################
library(lavaan)
library(tidySEM)
library(ggplot2)
library(vegan)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(semPlot)

#####################################
### Step 1: Load tables
#####################################

# ARG from 16S normalized table
ARG_data <- read_csv("Master/Metatranscriptomics/Excel_lists/ARG/ARG_data_combined.csv")

# Table RNA rarefied
data_rarefied_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_RNA.csv")

# Table Virus rarefied
data_rarefied_refseq <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_refseq_486.csv")


# Counttable (for Gernera identification)
count_data_pr2 <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_pr2.csv")
count_data_silva <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_silva.csv")

# Environmental tables
Environmental_Data_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metatranscriptomics.csv")
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")

# PCoA result Virus k=3
PCoA_Virus_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Virus_rarefied.csv")

# PCoA result Eukaryotes k=2
#PCoA_Euk_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Euk_rarefied.csv")
# PCoA result Eukaryotes Relativized!!!! k=2
PCoA_Euk_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Euk_rarefied_Rel.csv")

# PCoA result Prokaryotes k=2
#PCoA_Prok_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Prok_rarefied.csv")
# PCoA result Prokaryotes k=2 relativized!!!!!
PCoA_Prok_rarefied <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_Prok_rarefied_Rel.csv")

# PCoA result Eukaryotes k=2
PCoA_ARG <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/PCoA_ARG.csv")


# Set SRR_ID as column
data_rarefied_refseq <- data_rarefied_refseq %>%
  rename(SRR_ID = ...1)

# Set SRR_ID as column
data_rarefied_RNA <- data_rarefied_RNA %>%
  rename(SRR_ID = ...1)

#####################################
### Step 3: Remove outliers
#####################################

# Remove outlier samples
removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)

# Remove outlier samples
removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)
data_rarefied_refseq <- data_rarefied_refseq %>% filter(!SRR_ID %in% removed_refseq_samples)


# Remove outlier sample ARG DNA (missing preliminary samples)
#removed_ARG_samples <- c("SRR1611149", "SRR1544596", "SRR1611146")
#Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_ARG_samples)
#data_rarefied_refseq <-data_rarefied_refseq %>% filter(!SRR_ID %in% removed_ARG_samples)
#PCoA_Virus_rarefied <- PCoA_Virus_rarefied %>% filter(!SRR_ID %in% removed_ARG_samples)


# Remove outlier sample ARG RNA (missing preliminary samples)
#removed_ARG_samples_RNA <- c("SRR1611150", "SRR1544599", "SRR1611147")
#Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_ARG_samples_RNA)
#data_rarefied_RNA <- data_rarefied_RNA %>% filter(!SRR_ID %in% removed_ARG_samples_RNA)
#PCoA_Euk_rarefied <- PCoA_Euk_rarefied %>% filter(!SRR_ID %in% removed_ARG_samples_RNA)
#PCoA_Prok_rarefied <- PCoA_Prok_rarefied %>% filter(!SRR_ID %in% removed_ARG_samples_RNA) 


# select only neccessary Metadata columns
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten

Environmental_Data_DNA <- Environmental_Data_DNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten



#####################################
### Step 2: Prepare all tables and create Matrix
#####################################

#### Virus
# Set SRR_ID as rowname and convert to Matrix
virus_data_matrix <- data_rarefied_refseq %>%
  column_to_rownames(var = "SRR_ID") %>%
  as.matrix()
# relativize
virus_data_rel <- virus_data_matrix / rowSums(virus_data_matrix)

#### ARG
# Konvertiere ins Long-Format
ARG_long <- ARG_data %>%
  pivot_longer(cols = -Family, names_to = "SRR_ID", values_to = "Abundance")

# Konvertiere zurück ins Wide-Format, aber mit SRR_ID als Zeilenname
ARG_transformed <- ARG_long %>%
  pivot_wider(names_from = Family, values_from = Abundance, values_fill = 0)

removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
ARG_transformed <- ARG_transformed %>% filter(!SRR_ID %in% removed_refseq_samples)

# Set SRR_ID as rowname and convert to Matrix
ARG_matrix <- ARG_transformed %>%
  column_to_rownames(var = "SRR_ID") %>%
  as.matrix()

# relativize data to samples 
ARG_matrix <- ARG_matrix / rowSums(ARG_matrix)



#### Prokaryotes
# Extract unique Prokaryotic genus names from count_data_silva
prokaryote_genera <- unique(count_data_silva$Genus)

# Keep only prokaryotic genera that actually exist in data_rarefied_RNA
prokaryote_genera_present <- prokaryote_genera[prokaryote_genera %in% colnames(data_rarefied_RNA)]

# Now select only prokaryotic taxa from `data_rarefied_RNA`
data_prokaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(prokaryote_genera_present))  

# Set SRR_ID as rowname and convert to Matrix
Prok_data_matrix <- data_prokaryote %>%
  column_to_rownames(var = "SRR_ID") %>%
  as.matrix()

# relativize data to samples 
Prok_data_rel <- Prok_data_matrix / rowSums(Prok_data_matrix)



#### Eukaryotes
# Extract unique Eukaryotic genus names from count_data_pr2
eukaryote_genera <- unique(count_data_pr2$Genus)

# Keep only prokaryotic genera that actually exist in data_rarefied_RNA
eukaryote_genera_present <- eukaryote_genera[eukaryote_genera %in% colnames(data_rarefied_RNA)]

# Now select only prokaryotic taxa from `data_rarefied_RNA`
data_eukaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(eukaryote_genera_present))

# Set SRR_ID as rowname and convert to Matrix
euk_data_matrix <- data_eukaryote %>%
  column_to_rownames(var = "SRR_ID") %>%
  as.matrix()

# relativize data to samples 
euk_data_rel <- euk_data_matrix / rowSums(euk_data_matrix)

#####################################
### Step 3: Calculate Shannon Diversity
#####################################

#### Virus
Shannon_Vir <- diversity(virus_data_rel, index = "shannon")
# Shannon in eine neue Tabelle mit `SRR_ID` speichern
shannon_vir_table <- data.frame(SRR_ID = rownames(virus_data_rel), Shannon_Vir = Shannon_Vir)

#### ARG
Shannon_ARG <- diversity(ARG_matrix, index = "shannon")
shannon_arg_table <- data.frame(SRR_ID = rownames(ARG_matrix), Shannon_ARG = Shannon_ARG)

#### Prokaryotes
Shannon_Prok <- diversity(Prok_data_rel, index = "shannon")
shannon_prok_table <- data.frame(SRR_ID = rownames(Prok_data_rel), Shannon_Prok = Shannon_Prok)

#### Eukaryotes
Shannon_Euk <- diversity(euk_data_rel, index = "shannon")
shannon_euk_table <- data.frame(SRR_ID = rownames(euk_data_rel), Shannon_Euk = Shannon_Euk)


#####################################
### Step 4: Combine Shannon tables with Env Data
#####################################

#### Virus
shannon_vir_table <- shannon_vir_table %>%
  left_join(Environmental_Data_DNA %>% select(SRR_ID, Date), by = "SRR_ID") %>%
  select(-SRR_ID)

# combine PCoA result with Environmental table
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  left_join(shannon_vir_table, by = "Date")



#### ARG
shannon_arg_table <- shannon_arg_table %>%
  left_join(Environmental_Data_DNA %>% select(SRR_ID, Date), by = "SRR_ID") %>%
  select(-SRR_ID)

# combine PCoA result with Environmental table
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  left_join(shannon_arg_table, by = "Date")



#### Prokaryotes
# add Date to NMDS to combine NMDS and RNA rarefied data over Date and nor SRR_ID because this is different
shannon_prok_table <- shannon_prok_table %>%
  left_join(Environmental_Data_RNA %>% select(SRR_ID, Date), by = "SRR_ID") %>%
  select(-SRR_ID)

# combine PCoA result with Environmental table
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  left_join(shannon_prok_table, by = "Date")



#### Eukaryotes
# add Date to NMDS to combine NMDS and RNA rarefied data over Date and nor SRR_ID because this is different
shannon_euk_table <- shannon_euk_table %>%
  left_join(Environmental_Data_RNA %>% select(SRR_ID, Date), by = "SRR_ID") %>%
  select(-SRR_ID)

# combine PCoA result with Environmental table
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  left_join(shannon_euk_table, by = "Date")



#####################################
### Step 4: Combine PCoA Results with Env Data
#####################################
# add Date to NMDS to combine NMDS and RNA rarefied data over Date and not SRR_ID because this is different
PCoA_Virus_rarefied <- PCoA_Virus_rarefied %>%
  left_join(Environmental_Data_DNA %>% select(SRR_ID, Date), by = "SRR_ID")

# rename PC1 coumn that it is clear it comes from Virus
PCoA_Virus_rarefied <- PCoA_Virus_rarefied %>%
  select(Date, PC1, PC2) %>%
  rename(Vir_PC1= PC1, Vir_PC2 = PC2)

# combine PCoA result with Environmental table
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  left_join(PCoA_Virus_rarefied, by = "Date")



# add Date to NMDS to combine NMDS and RNA rarefied data over Date and not SRR_ID because this is different
PCoA_ARG <- PCoA_ARG %>%
  left_join(Environmental_Data_DNA %>% select(SRR_ID, Date), by = "SRR_ID")

# rename PC1 coumn that it is clear it comes from Virus
PCoA_ARG <- PCoA_ARG %>%
  select(Date, PC1, PC2) %>%
  rename(ARG_PC1= PC1, ARG_PC2 = PC2)

# combine PCoA result with Environmental table
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  left_join(PCoA_ARG, by = "Date")



# add Date to NMDS to combine NMDS and RNA rarefied data over Date and nor SRR_ID because this is different
PCoA_Euk_rarefied <- PCoA_Euk_rarefied %>%
  left_join(Environmental_Data_RNA %>% select(SRR_ID, Date), by = "SRR_ID")

# rename PC1 coumn that it is clear it comes from Virus
PCoA_Euk_rarefied <- PCoA_Euk_rarefied %>%
  select(Date, PC1, PC2) %>%
  rename(Euk_PC1 = PC1, Euk_PC2 = PC2)

# combine PCoA result with Environmental table
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  left_join(PCoA_Euk_rarefied, by = "Date")



# add Date to NMDS to combine NMDS and RNA rarefied data over Date and nor SRR_ID because this is different
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
### Step 4: Create one table with all data
#####################################
# Dataframe mit allen relevanten Variablen für SEM erstellen
sem_data <- Environmental_Data_DNA




#####################################
### Step 5: Standardize shannon with Scale 
#####################################
# To make the Shannon Values from Prok, Euk, ARG, and Virus comparable! Because they are very different.
#sem_data <- sem_data %>%
#mutate(
#Shannon_Vir = as.numeric(scale(Shannon_Vir)),
#Shannon_ARG = as.numeric(scale(Shannon_ARG)),
#Shannon_Prok = as.numeric(scale(Shannon_Prok)),
#Shannon_Euk = as.numeric(scale(Shannon_Euk))
#)


#####################################
### Step 6: Construct SEM Model
#####################################
sem_data <- sem_data %>%
  mutate(Season = as.numeric(factor(Season, level=c( "spring", "summer", "autumn", "winter"))))  # Faktor in Zahlen umwandeln


# Test Model mit Shannon für alle
sem_model <- '
  # Einfluss der Season auf Euk, Prok, Vir
  Shannon_Euk ~ Season
  Shannon_Vir ~ Season
  Shannon_Prok ~ Season
  
  # Einfluss der Viren auf Prokaryoten
  Shannon_Prok ~ Shannon_Vir + Shannon_Euk
  
    # Einfluss der Euk auf Viren
  Shannon_Euk ~ Shannon_Vir
  
  # Einfluss der Viren auf ARG
  Shannon_ARG ~ Shannon_Vir + Shannon_Prok + Shannon_Euk
  
'


# SEM Modell fitten
fit <- sem(sem_model, data = sem_data)


# Ergebnisse ausgeben
summary(fit, fit.measures = TRUE, rsquare = TRUE, standardized = TRUE)

# Define Layout
sem_layout <- get_layout(
  NA , "Season", NA,
  "Shannon_Vir", NA, "Shannon_Euk", 
  "Shannon_ARG" , NA, "Shannon_Prok",
  rows = 3
)



##########################
### Model with Estimate
##########################
# Öffne PDF-Gerät
pdf("sem_plot_shannon_estimate.pdf", width = 10, height = 8)

# Zeichne den Plot mit unstandardisierten Estimates
semPaths(
  fit,
  whatLabels = "est",  # <- hier ist der Unterschied!
  layout = "tree",
  edge.label.cex = 0.8,
  sizeMan = 8,
  sizeLat = 10,
  nCharNodes = 0,
  curvePivot = TRUE
)

# Schließe das PDF-Gerät
dev.off()


