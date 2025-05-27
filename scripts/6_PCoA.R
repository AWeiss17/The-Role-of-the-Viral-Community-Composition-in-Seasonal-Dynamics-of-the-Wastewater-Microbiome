#################################################################################################
########################################### Viral PCoA ##########################################
#################################################################################################


####################################
### Load Required Libraries
#####################################

library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)
library(readr)
library(ggalt)


#####################################
### Step 1: Load rarefied counttables & environmental tables
#####################################

data_rarefied_refseq <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_refseq_486.csv")

Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")

# Set SRR_ID as column
data_rarefied_refseq <- data_rarefied_refseq %>%
  rename(SRR_ID = ...1)

Environmental_Data_DNA <- Environmental_Data_DNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten

# Remove outlier samples
removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
data_rarefied_refseq <- data_rarefied_refseq %>% filter(!SRR_ID %in% removed_refseq_samples)
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)

# Liste der Virusfamilien, die entfernt werden sollen
#families_to_remove <- c("Iridoviridae", "Mimiviridae")

# Entferne diese Familien
#data_rarefied_refseq <- data_rarefied_refseq %>%
#select(-one_of(families_to_remove))


#####################################
### Step 2: Set up distance Matrix for PCoA calculation
#####################################
# PCoA needs distance matrix like Bray-Curtis

# Set SRR_ID as rowname and convert to Matrix
virus_data_matrix <- data_rarefied_refseq %>%
  column_to_rownames(var = "SRR_ID") %>%
  as.matrix()
# Normalization of the data is not neccessary because we use the rarefied table, which is already "normalized"
# Calculate Bray-Curtis-Disssimilarity Matrix
bray_curtis <- vegdist(virus_data_matrix, method = "bray")


#####################################
### Step 3: Calculate PCoA
#####################################
# Conduct PCoA
# cmdscale = shows samples based on distances (unconstrained of Metadata)
# good for only visualizing dissimilarities / similarities bewteen samples
pcoa_result <- cmdscale(bray_curtis, eig = TRUE, k = 2)  # k = 2 für 2D-Darstellung

# Set up DataFrame with PCoA-Coordinates
pcoa_df <- data.frame(SRR_ID = rownames(pcoa_result$points),
                      PC1 = pcoa_result$points[, 1],
                      PC2 = pcoa_result$points[, 2])


write_csv(pcoa_df, "PCoA_Virus_new.csv")
# Calculate Variance of axis
eig_values <- pcoa_result$eig
variance_explained <- eig_values / sum(eig_values) * 100
cat("Erklärte Varianz der ersten beiden Achsen:\n")
cat("PC1:", round(variance_explained[1], 2), "%\n")
cat("PC2:", round(variance_explained[2], 2), "%\n")

# Scree plot
# Shows Eigenwert-Diagram, and show how much variance is explained by each axis (PC1, PC2, PC3,...PC45)
# Visualization shows: First two axes (PC1 & PC2) explain biggest part of variation in data
# this says that 2D visualisation is sufficient
barplot(eig_values, names.arg = 1:length(eig_values), 
        main = "Scree Plot", xlab = "Hauptachsen", ylab = "Eigenwerte")

#####################################
### Step 4: Add Environmental Data
#####################################

# Add Environmental Data to PCoA Coordinates
pcoa_env_data <- pcoa_df %>%
  left_join(Environmental_Data_DNA, by = "SRR_ID")



#####################################
### PCoA over Seasons
#####################################

# Define correct order of seasons in legend
pcoa_env_data$Season <- factor(pcoa_env_data$Season, levels = c("spring", "summer", "autumn", "winter"))

# Farben für Seasons
season_colors <- c("spring" = "#62b6cb",
                   "summer" = "#ff7d00",
                   "autumn" = "#8c1c13", 
                   "winter" = "#134074")
# PCoA plot mit ggplot2 (färbt nach Season)
ggplot(pcoa_env_data, aes(x = PC1, y = PC2, color = Season)) +
  geom_point(size = 3) +
  labs(title = "PCoA - Bray-Curtis Dissimilarity",
       x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 2), "%)")) +
  theme_minimal() +
  scale_color_manual(values = season_colors)


ggsave(
  filename = "PCoA_Virus_Season_rarefied.png",  # Output file name
  plot = last_plot(),  # Save the last generated plot
  device = "png",  # Save as PNG
  width = 10,  # Adjust width (in inches)
  height = 7,  # Adjust height (in inches)
  units = "in",  # Units for width and height
  dpi = 600,  # High resolution (300-600 DPI recommended for publication)
  bg = "white" # Set backgrund white
)


## plot with polygons
# PCoA plot mit ggplot2 (färbt nach Season)
ggplot(pcoa_env_data, aes(x = PC1, y = PC2, color = Season)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_encircle(aes(fill = Season), alpha = 0.2, s_shape = 1, expand = 0.02) +  # Polygone für Gruppen
  labs(title = "PCoA ",
       x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 2), "%)")) +
  theme_minimal() +
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +  # Gleiche Farben für die Füllung
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )

ggsave(
  filename = "PCoA_Vir.png",  
  plot = last_plot(),  
  device = "png",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)

library(svglite)
ggsave(
  filename = "PCoA_Vir.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)


######################################################################################################
########################################## Prokaryotes PCoA ##########################################
######################################################################################################


####################################
### Load Required Libraries
#####################################

library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)
library(readr)
library(ggalt)


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


#####################################
### Step 2: Extract Prokaryotes & Eukaryotes from rarefied table
#####################################

# Set SRR_ID as column
data_rarefied_RNA <- data_rarefied_RNA %>%
  rename(SRR_ID = ...1)


# Extract Prokaryotes
# Extract unique Prokaryotic genus names from count_data_silva
prokaryote_genera <- unique(count_data_silva$Genus)

# Keep only prokaryotic genera that actually exist in data_rarefied_RNA
prokaryote_genera_present <- prokaryote_genera[prokaryote_genera %in% colnames(data_rarefied_RNA)]

# Now select only prokaryotic taxa from `data_rarefied_RNA`
data_prokaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(prokaryote_genera_present))  


# Extract Eukaryotes
# Extract unique Eukaryotic genus names from count_data_pr2
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

# select only neccessary Metadata columns
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten


#####################################
### Step 2: Set up distance Matrix for PCoA calculation
#####################################
# PCoA needs distance matrix like Bray-Curtis

# Set SRR_ID as rowname and convert to Matrix
Prok_data_matrix <- data_prokaryote %>%
  column_to_rownames(var = "SRR_ID") %>%
  as.matrix()

#####################################
### Step 6: Hellinger Transformation
#####################################
# helps to mimimalize rare taxa effect 
#data_prokaryote_hel <- decostand(Prok_data_matrix, method = "hellinger")
#bray_curtis <- vegdist(data_prokaryote_hel, method = "bray")


# Normalization of the data is neccessary!
# Although we use the rarefied table, which is already "normalized" 
# BUT! We divide pr2 and silva dataset! So new relativization is neccessary!
# Relativieren pro Probe (Zeile)
Prok_data_rel <- Prok_data_matrix / rowSums(Prok_data_matrix)


# Calculate Bray-Curtis-Disssimilarity Matrix
bray_curtis <- vegdist(Prok_data_rel, method = "bray")


#####################################
### Step 3: Calculate PCoA
#####################################
# Conduct PCoA
# cmdscale = shows samples based on distances (unconstrained of Metadata)
# good for only visualizing dissimilarities / similarities bewteen samples
pcoa_result <- cmdscale(bray_curtis, eig = TRUE, k = 2)  # k = 2 für 2D-Darstellung


# Set up DataFrame with PCoA-Coordinates
pcoa_df <- data.frame(SRR_ID = rownames(pcoa_result$points),
                      PC1 = pcoa_result$points[, 1],
                      PC2 = pcoa_result$points[, 2])

# Calculate Variance of axis
eig_values <- pcoa_result$eig
variance_explained <- eig_values / sum(eig_values) * 100
cat("Erklärte Varianz der ersten beiden Achsen:\n")
cat("PC1:", round(variance_explained[1], 2), "%\n")
cat("PC2:", round(variance_explained[2], 2), "%\n")

# Scree plot
# Shows Eigenwert-Diagram, and show how much variance is explained by each axis (PC1, PC2, PC3,...PC45)
# Visualization shows: First two axes (PC1 & PC2) explain biggest part of variation in data
# this says that 2D visualisation is sufficient
barplot(eig_values, names.arg = 1:length(eig_values), 
        main = "Scree Plot", xlab = "Hauptachsen", ylab = "Eigenwerte")



#####################################
### Step 4: Add Environmental Data
#####################################

# Add Environmental Data to PCoA Coordinates
pcoa_env_data <- pcoa_df %>%
  left_join(Environmental_Data_RNA, by = "SRR_ID")



#####################################
### PCoA over Seasons
#####################################

# Define correct order of seasons in legend
pcoa_env_data$Season <- factor(pcoa_env_data$Season, levels = c("spring", "summer", "autumn", "winter"))


season_colors <- c("spring" = "#62b6cb",
                   "summer" = "#ff7d00",
                   "autumn" = "#8c1c13", 
                   "winter" = "#134074")
# PCoA plot mit ggplot2 (färbt nach Season)
ggplot(pcoa_env_data, aes(x = PC1, y = PC2, color = Season)) +
  geom_point(size = 3) +
  labs(title = "PCoA",
       x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 2), "%)")) +
  theme_minimal() +
  scale_color_manual(values = season_colors)


## plot with polygons
# PCoA plot mit ggplot2 (färbt nach Season)
ggplot(pcoa_env_data, aes(x = PC1, y = PC2, color = Season)) +
  geom_point(size = 3) +
  geom_encircle(aes(fill = Season), alpha = 0.2, s_shape = 1, expand = 0.02) +  
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  labs(title = "PCoA",
       x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 2), "%)")) +
  theme_minimal() +
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )

ggsave(
  filename = "PCoA_Prok.png",  
  plot = last_plot(),  
  device = "png",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)

library(svglite)
ggsave(
  filename = "PCoA_Prok.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)




#####################################################################################################
########################################## Eukaryotes PCoA ##########################################
#####################################################################################################


####################################
### Load Required Libraries
#####################################

library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)
library(readr)
library(ggalt)


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


#####################################
### Step 2: Extract Prokaryotes & Eukaryotes from rarefied table
#####################################

# Set SRR_ID as column
data_rarefied_RNA <- data_rarefied_RNA %>%
  rename(SRR_ID = ...1)


# Extract Prokaryotes
# Extract unique Prokaryotic genus names from count_data_silva
prokaryote_genera <- unique(count_data_silva$Genus)

# Keep only prokaryotic genera that actually exist in data_rarefied_RNA
prokaryote_genera_present <- prokaryote_genera[prokaryote_genera %in% colnames(data_rarefied_RNA)]

# Now select only prokaryotic taxa from `data_rarefied_RNA`
data_prokaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(prokaryote_genera_present))  


# Extract Eukaryotes
# Extract unique Eukaryotic genus names from count_data_pr2
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

# select only neccessary Metadata columns
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten


#####################################
### Step 2: Set up distance Matrix for PCoA calculation
#####################################
# PCoA needs distance matrix like Bray-Curtis

# Set SRR_ID as rowname and convert to Matrix
euk_data_matrix <- data_eukaryote %>%
  column_to_rownames(var = "SRR_ID") %>%
  as.matrix()

# Normalization of the data is neccessary!
# Although we use the rarefied table, which is already "normalized" 
# BUT! We divide pr2 and silva dataset! So new relativization is neccessary!
# Relativieren pro Probe (Zeile)
euk_data_rel <- euk_data_matrix / rowSums(euk_data_matrix)


# Calculate Bray-Curtis-Disssimilarity Matrix
bray_curtis <- vegdist(euk_data_rel, method = "bray")


#####################################
### Step 3: Calculate PCoA
#####################################
# Conduct PCoA
# cmdscale = shows samples based on distances (unconstrained of Metadata)
# good for only visualizing dissimilarities / similarities bewteen samples
pcoa_result <- cmdscale(bray_curtis, eig = TRUE, k = 2)  # k = 2 für 2D-Darstellung

# Set up DataFrame with PCoA-Coordinates
pcoa_df <- data.frame(SRR_ID = rownames(pcoa_result$points),
                      PC1 = pcoa_result$points[, 1],
                      PC2 = pcoa_result$points[, 2])

# Calculate Variance of axis
eig_values <- pcoa_result$eig
variance_explained <- eig_values / sum(eig_values) * 100
cat("Erklärte Varianz der ersten beiden Achsen:\n")
cat("PC1:", round(variance_explained[1], 2), "%\n")
cat("PC2:", round(variance_explained[2], 2), "%\n")

# Scree plot
# Shows Eigenwert-Diagram, and show how much variance is explained by each axis (PC1, PC2, PC3,...PC45)
# Visualization shows: First two axes (PC1 & PC2) explain biggest part of variation in data
# this says that 2D visualisation is sufficient
barplot(eig_values, names.arg = 1:length(eig_values), 
        main = "Scree Plot", xlab = "Hauptachsen", ylab = "Eigenwerte")



#####################################
### Step 4: Add Environmental Data
#####################################

# Add Environmental Data to PCoA Coordinates
pcoa_env_data <- pcoa_df %>%
  left_join(Environmental_Data_RNA, by = "SRR_ID")



#####################################
### PCoA over Seasons
#####################################
season_colors <- c("spring" = "#62b6cb",
                   "summer" = "#ff7d00",
                   "autumn" = "#8c1c13", 
                   "winter" = "#134074")
# Define correct order of seasons in legend
pcoa_env_data$Season <- factor(pcoa_env_data$Season, levels = c("spring", "summer", "autumn", "winter"))

# PCoA plot mit ggplot2 (färbt nach Season)
ggplot(pcoa_env_data, aes(x = PC1, y = PC2, color = Season)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  labs(title = "PCoA",
       x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 2), "%)")) +
  theme_minimal() +
  scale_color_manual(values = season_colors) +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )

## plot with polygons
# PCoA plot mit ggplot2 (färbt nach Season)
ggplot(pcoa_env_data, aes(x = PC1, y = PC2, color = Season)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_encircle(aes(fill = Season), alpha = 0.2, s_shape = 1, expand = 0.02) +  # Polygone für Gruppen
  labs(title = "PCoA",
       x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 2), "%)")) +
  theme_minimal() +
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )


ggsave(
  filename = "PCoA_Euk.png",  
  plot = last_plot(),  
  device = "png",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)

library(svglite)
ggsave(
  filename = "PCoA_Euk.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)

###############################################################################################
########################################### ARG PCoA ##########################################
###############################################################################################


####################################
### Load Required Libraries
#####################################

library(vegan)
library(ggplot2)
library(dplyr)
library(tibble)
library(tidyr)
library(readr)
library(ggalt)
library(ape)

#####################################
### Step 1: Load rarefied counttables & environmental tables
#####################################

# ARG from 16S normalized table
#ARG_data <- read_csv("Master/Metatranscriptomics/Excel_lists/ARG/Family_normalized_16S.csv")
ARG_data <- read_csv("Master/Metatranscriptomics/Excel_lists/ARG/ARG_data_combined.csv")
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")
data_rarefied_refseq <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_refseq_486.csv")


Environmental_Data_DNA <- Environmental_Data_DNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten

# Set SRR_ID as column
data_rarefied_refseq <- data_rarefied_refseq %>%
  rename(SRR_ID = ...1)


#####################################
### Step 2: Prepare ARG table
#####################################

# Konvertiere ins Long-Format
ARG_long <- ARG_data %>%
  pivot_longer(cols = -Family, names_to = "SRR_ID", values_to = "Abundance")

# Konvertiere zurück ins Wide-Format, aber mit SRR_ID als Zeilenname
ARG_transformed <- ARG_long %>%
  pivot_wider(names_from = Family, values_from = Abundance, values_fill = 0)


#####################################
### Step 3: Remove outliers
#####################################

removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
ARG_transformed <- ARG_transformed %>% filter(!SRR_ID %in% removed_refseq_samples)
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)
data_rarefied_refseq <- data_rarefied_refseq %>% filter(!SRR_ID %in% removed_refseq_samples)

#removed_ARG_samples <- c("SRR1611149", "SRR1544596", "SRR1611146")
#Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_ARG_samples)

# Setze SRR_ID als Zeilenname für die NMDS-Analyse
data_rarefied_refseq <- data_rarefied_refseq %>% column_to_rownames(var = "SRR_ID")


#####################################
### Step 2: Set up distance Matrix for PCoA calculation
#####################################
# PCoA needs distance matrix like Bray-Curtis

# Set SRR_ID as rowname and convert to Matrix
ARG_matrix <- ARG_transformed %>%
  column_to_rownames(var = "SRR_ID") %>%
  as.matrix()

# relativize data to samples 
ARG_matrix <- ARG_matrix / rowSums(ARG_matrix)

# Calculate Bray-Curtis-Disssimilarity Matrix
bray_curtis <- vegdist(ARG_matrix, method = "bray")


#####################################
### Step 3: Calculate PCoA
#####################################
# Conduct PCoA
# cmdscale = shows samples based on distances (unconstrained of Metadata)
# good for only visualizing dissimilarities / similarities bewteen samples
pcoa_result <- cmdscale(bray_curtis, eig = TRUE, k = 2)  # k = 2 für 2D-Darstellung

# Set up DataFrame with PCoA-Coordinates
pcoa_df <- data.frame(SRR_ID = rownames(pcoa_result$points),
                      PC1 = pcoa_result$points[, 1],
                      PC2 = pcoa_result$points[, 2])

# Calculate Variance of axis
eig_values <- pcoa_result$eig
variance_explained <- eig_values / sum(eig_values) * 100
cat("Erklärte Varianz der ersten beiden Achsen:\n")
cat("PC1:", round(variance_explained[1], 2), "%\n")
cat("PC2:", round(variance_explained[2], 2), "%\n")

# Scree plot
# Shows Eigenwert-Diagram, and show how much variance is explained by each axis (PC1, PC2, PC3,...PC45)
# Visualization shows: First two axes (PC1 & PC2) explain biggest part of variation in data
# this says that 2D visualisation is sufficient
barplot(eig_values, names.arg = 1:length(eig_values), 
        main = "Scree Plot", xlab = "Hauptachsen", ylab = "Eigenwerte")

#####################################
### Step 4: Add Environmental Data
#####################################

# Add Environmental Data to PCoA Coordinates
pcoa_env_data <- pcoa_df %>%
  left_join(Environmental_Data_DNA, by = "SRR_ID")


write_csv(pcoa_df, "PCoA_ARG.csv")
#####################################
### PCoA over Seasons
#####################################

# Define correct order of seasons in legend
pcoa_env_data$Season <- factor(pcoa_env_data$Season, levels = c("spring", "summer", "autumn", "winter"))

# Farben für Seasons
season_colors <- c("spring" = "#62b6cb",
                   "summer" = "#ff7d00",
                   "autumn" = "#8c1c13", 
                   "winter" = "#134074")


# PCoA plot
ggplot(pcoa_env_data, aes(x = PC1, y = PC2, color = Season)) +
  geom_point(size = 3) +
  geom_encircle(aes(fill = Season), alpha = 0.2, s_shape = 1, expand = 0.02) + 
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  labs(title = "PCoA",
       x = paste0("PC1 (", round(variance_explained[1], 2), "%)"),
       y = paste0("PC2 (", round(variance_explained[2], 2), "%)")) +
  theme_minimal() +
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +  
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )

ggsave(
  filename = "PCoA_ARG.png",  
  plot = last_plot(),  
  device = "png",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)

library(svglite)
ggsave(
  filename = "PCoA_ARG.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)
