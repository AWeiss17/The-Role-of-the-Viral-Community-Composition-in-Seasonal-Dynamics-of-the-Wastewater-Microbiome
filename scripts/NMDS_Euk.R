#####################################################################################################
########################################## Eukaryota NMDS ##########################################
#####################################################################################################


####################################
### Load Required Libraries
#####################################

library(readr)
library(ggalt)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)



# Table RNA rarefied
data_rarefied_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_RNA.csv")

# Counttable (for Gernera identification)
count_data_pr2 <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_pr2.csv")
count_data_silva <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_silva.csv")

# Environmental RNA table
Environmental_Data_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metatranscriptomics.csv")
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")


#####################################
### Step 2: Extract Eukaryotes from rarefied table
#####################################

# Set SRR_ID as column
data_rarefied_RNA <- data_rarefied_RNA %>%
  rename(SRR_ID = ...1)


#####################################
### Step 3: Extract Eukaryotes Genera & map at Order
#####################################
# Nur Genera auswählen, die in data_rarefied_RNA vorkommen
eukaryote_genera_present <- unique(count_data_pr2$Genus[count_data_pr2$Genus %in% colnames(data_rarefied_RNA)])

# Filtere count_data_pr2, um nur relevante Ordnungen zu behalten
genus_to_order_pr2_filtered <- count_data_pr2 %>%
  filter(Genus %in% eukaryote_genera_present) %>%
  select(Genus, Order) %>%
  distinct()

# Daten der relevanten Genera aus der RNA-Tabelle extrahieren
data_eukaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(eukaryote_genera_present))

# Genera-Spalten zu ihrer zugehörigen Ordnung umbenennen, aber doppelte vermeiden
genus_to_order_map <- setNames(genus_to_order_pr2_filtered$Order, genus_to_order_pr2_filtered$Genus)

# ️Spalten korrekt umbenennen, indem wir eine Hilfstabelle erstellen
data_eukaryote_order <- data_eukaryote %>%
  pivot_longer(-SRR_ID, names_to = "Genus", values_to = "Abundance") %>%
  left_join(genus_to_order_pr2_filtered, by = "Genus") %>%
  group_by(SRR_ID, Order) %>%  # Gruppiere nach SRR_ID und Order (nicht Genus!)
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Order, values_from = Abundance, values_fill = 0)  # Zurück ins `wide`-Format



#####################################
### Step 3: Remove Outlier from Metadata
#####################################

# Remove outlier samples (are already exluded from rarefied data table, son only exclude them from Env. table)
removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)

# select only neccessary Metadata columns
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten



# Setze SRR_ID als Zeilenname für die NMDS-Analyse
data_eukaryote_order <- data_eukaryote_order %>% column_to_rownames(var = "SRR_ID")



#####################################
### Step 2: Choose optimal number of dimensions
#####################################
# Normalization of the data is neccessary!
# Although we use the rarefied table, which is already "normalized" 
# BUT! We divide pr2 and silva dataset! So new relativization is neccessary!
# Relativieren pro Probe (Zeile)
euk_data_rel <- data_eukaryote_order / rowSums(data_eukaryote_order)
# Perform NMDS with 1 to 4 dimensions and save the stress value
stress <- c()
for (i in 1:4) {
  stress <- c(stress, metaMDS(euk_data_rel, dist="bray", k = i)$stress)
}

# Plot the stress values vs dimensions   
plot(1:4, stress, xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")

# Save stress plot
png("Euk_NMDS_StressPlot.png", width = 800, height = 600, res = 150)
plot(1:4, stress, xlab = "# of Dimensions", ylab = "Stress", main = "Virus NMDS stress plot")
dev.off()

#####################################
### Step 3: Calculate NMDS
#####################################
# As the final result depends on the initial random placement of the points 
# we set a seed to make the results reproducible
set.seed(183672807)

# Carry out NMDS
NMDS <- metaMDS(euk_data_rel, dist="bray", k=2, maxit=1000, trymax=100, autotransform = FALSE)
# k = number of dimensions (in this case I choose 2 Dimensions, because it refers to a stress value of a good fit ~0.15)
# trymax = maximum numbers of random starts in search of stable solution
# maxit = Maximum number of iterations in the single NMDS run
# autotransform = applies a square root transformation on the community matrix

# Plot NMDS
plot(NMDS, display = "sites")

Stress <- NMDS$stress
head(Stress)

# NMDS-Koordinaten extrahieren
NMDS_Euk_df <- as.data.frame(scores(NMDS, display = "sites")) %>%
  rownames_to_column("SRR_ID")  # SRR_ID als Spalte hinzufügen

# Speichern der NMDS-Koordinaten als CSV
write_csv(NMDS_Euk_df, "NMDS_Euk_Coordinates.csv")
#####################################
### Step 5: Fit environmental variables onto the NMDS of refseq
#####################################

# Convert tibble to data frame if needed
Environmental_Data_RNA <- as.data.frame(Environmental_Data_RNA)

# Set `SRR_ID` as row name
rownames(Environmental_Data_RNA) <- Environmental_Data_RNA$SRR_ID

# Remove the `SRR_ID` column
Environmental_Data_RNA <- Environmental_Data_RNA[, -1]

# Stelle sicher, dass Reihenfolge der SRR_IDs in beiden Datensätzen gleich ist
Environmental_Data_RNA <- Environmental_Data_RNA[match(rownames(euk_data_rel), rownames(Environmental_Data_RNA)), ]

# Überprüfe, ob Reihenfolge übereinstimmt
identical(rownames(euk_data_rel), rownames(Environmental_Data_RNA)) 

# Check if there are any NA values in the dataset
anyNA(Environmental_Data_RNA)

# Fit environmental vectors and factors onto the ordination
env <- envfit(NMDS, Environmental_Data_RNA, permutations = 999)

# Plot environmental vectors
plot(NMDS, display = "sites")
plot(env, p.max = 0.05)


#####################################
### Step 6: Extract only significant environmental facors
#####################################

# Extract p-values from environmental factors
env_pvalues <- env$vectors$pvals

# Filter only significant Vectors (p < 0.05)
significant_env_vectors <- as.data.frame(scores(env, "vectors")) * ordiArrowMul(env)
significant_env_vectors$Feature <- rownames(significant_env_vectors)
significant_env_vectors <- significant_env_vectors[env_pvalues < 0.05, ]

# Erstelle eine neue Spalte mit den gewünschten Bezeichnungen
label_mapping <- c("pH_online" = "pH", "NO3_online_mg.L" = "Nitrate", "Temp_manual" = "Water temperature", "Temp_Air" = "Air temperature")

# Neue Labels in die `significant_env_vectors`-Tabelle übernehmen
significant_env_vectors$Feature_renamed <- label_mapping[significant_env_vectors$Feature]

# Falls ein Feature nicht in `label_mapping` ist, behalte den Originalnamen
significant_env_vectors$Feature_renamed[is.na(significant_env_vectors$Feature_renamed)] <- significant_env_vectors$Feature


# Extract NMDS-Coordinates
nmds_points <- as.data.frame(scores(NMDS, display = "sites"))
nmds_points$Sample <- rownames(nmds_points)

# Add seasons
nmds_points$Season <- Environmental_Data_RNA$Season[match(rownames(nmds_points), rownames(Environmental_Data_RNA))]

# Season as factor
nmds_points$Season <- as.factor(nmds_points$Season)

# set order within legend
nmds_points$Season <- factor(nmds_points$Season, levels = c("spring", "summer", "autumn", "winter"))

# Define Season colors
season_colors <- c("spring" = "#62b6cb", 
                   "summer" = "#ff7d00", 
                   "autumn" = "#8c1c13", 
                   "winter" = "#134074")



# Extrahiere Umweltvektoren aus envfit
env_vectors <- as.data.frame(scores(env, "vectors")) * ordiArrowMul(env)
env_vectors$Feature <- rownames(env_vectors)

# Plot with Environmental vectors
ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Season)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_polygon(data = nmds_points %>% group_by(Season) %>% slice(chull(NMDS1, NMDS2)), 
               aes(fill = Season, group = Season), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_segment(data = significant_env_vectors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  geom_text_repel(data = significant_env_vectors, aes(x = NMDS1, y = NMDS2, label = Feature_renamed),
                  size = 4, color = "black") +
  annotate("text", x = Inf, y = -Inf, label = paste("Stress =", round(NMDS$stress, 3)), 
           hjust = 1, vjust = -0.5, size = 3.5) +
  
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  theme_minimal() +
  labs(title = "Eukaryota", x = "NMDS1", y = "NMDS2") +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )


ggsave(
  filename = "NMDS_Euk_Env.png",  
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
  filename = "NMDS_Euk_Env.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)

#####################################################################################################
########################################## Prokaryota NMDS with top 15 Eukaryotes as Vector ##########################################
#####################################################################################################


####################################
### Load Required Libraries
#####################################

library(readr)
library(ggalt)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(vegan)
library(dplyr)
library(ggrepel)
library(tidyr)
library(tibble)



# Table RNA rarefied
data_rarefied_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_RNA.csv")

# Counttable (for Gernera identification)
count_data_pr2 <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_pr2.csv")
count_data_silva <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_silva.csv")

# Environmental RNA table
Environmental_Data_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metatranscriptomics.csv")
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")


#####################################
### Step 2: Extract Prokaryotes & Eukaryotes from rarefied table
#####################################

# Set SRR_ID as column
data_rarefied_RNA <- data_rarefied_RNA %>%
  rename(SRR_ID = ...1)


#####################################
### Step 3: Extract Eukaryotes Genera & map at Order
#####################################
# Nur Genera auswählen, die in data_rarefied_RNA vorkommen
eukaryote_genera_present <- unique(count_data_pr2$Genus[count_data_pr2$Genus %in% colnames(data_rarefied_RNA)])

# Filtere count_data_pr2, um nur relevante Ordnungen zu behalten
genus_to_order_pr2_filtered <- count_data_pr2 %>%
  filter(Genus %in% eukaryote_genera_present) %>%
  select(Genus, Order) %>%
  distinct()

# Daten der relevanten Genera aus der RNA-Tabelle extrahieren
data_eukaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(eukaryote_genera_present))

# Genera-Spalten zu ihrer zugehörigen Ordnung umbenennen, aber doppelte vermeiden
genus_to_order_map <- setNames(genus_to_order_pr2_filtered$Order, genus_to_order_pr2_filtered$Genus)

# ️Spalten korrekt umbenennen, indem wir eine Hilfstabelle erstellen
data_eukaryote_order <- data_eukaryote %>%
  pivot_longer(-SRR_ID, names_to = "Genus", values_to = "Abundance") %>%
  left_join(genus_to_order_pr2_filtered, by = "Genus") %>%
  group_by(SRR_ID, Order) %>%  # Gruppiere nach SRR_ID und Order (nicht Genus!)
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Order, values_from = Abundance, values_fill = 0)  # Zurück ins `wide`-Format



#####################################
### Step 3: Remove Outlier from Metadata
#####################################

# Remove outlier samples (are already exluded from rarefied data table, son only exclude them from Env. table)
removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)

# Select only necessary Metadata columns
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season)  

# Set SRR_ID as rowname
data_eukaryote_order <- data_eukaryote_order %>% column_to_rownames(var = "SRR_ID")


#####################################
### Step 2: Choose optimal number of dimensions
#####################################
# Normalization of the data is neccessary!
# Although we use the rarefied table, which is already "normalized" 
# BUT! We divide pr2 and silva dataset! So new relativization is neccessary!
# Relativieren pro Probe (Zeile)
euk_data_rel <- data_eukaryote_order / rowSums(data_eukaryote_order)
# Perform NMDS with 1 to 4 dimensions and save the stress value
stress <- c()
for (i in 1:4) {
  stress <- c(stress, metaMDS(euk_data_rel, dist="bray", k = i)$stress)
}

# Plot the stress values vs dimensions   
plot(1:4, stress, xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")


#####################################
### Step 3: Calculate NMDS
#####################################
# As the final result depends on the initial random placement of the points 
# we set a seed to make the results reproducible
set.seed(183672807)

# Carry out NMDS
NMDS <- metaMDS(euk_data_rel, dist="bray", k=2, maxit=1000, trymax=100, autotransform = FALSE)
# k = number of dimensions (in this case I choose 2 Dimensions, because it refers to a stress value of a good fit ~0.15)
# trymax = maximum numbers of random starts in search of stable solution
# maxit = Maximum number of iterations in the single NMDS run
# autotransform = applies a square root transformation on the community matrix

# Plot NMDS
plot(NMDS, display = "sites")

Stress <- NMDS$stress
head(Stress)


#####################################
### Step 6: Fit Prokaryote-Taxa as Vectors
#####################################

# Extract Top 15 abundant Taxa
top_eukaryotes <- colSums(euk_data_rel) %>%
  sort(decreasing = TRUE) %>%
  head(15)

# Filter data on top 15
eukaryote_data_top <- euk_data_rel[, names(top_eukaryotes)]

# Fit NMDS with Prokaryotes-Taxa as Vectors
eukaryote_vectors <- envfit(NMDS, eukaryote_data_top, permutations = 999)

# Extract p-values and significant Taxa (p<0.05)
eukaryote_pvalues <- eukaryote_vectors$vectors$pvals

# Behalte nur signifikante Taxa (p < 0.05)
significant_eukaryote_vectors <- as.data.frame(scores(eukaryote_vectors, "vectors")) * ordiArrowMul(eukaryote_vectors)
significant_eukaryote_vectors$Feature <- rownames(significant_eukaryote_vectors)
significant_eukaryote_vectors <- significant_eukaryote_vectors[eukaryote_pvalues < 0.05, ]


#####################################
### Step 7: NMDS-Plot erstellen
#####################################

# Extract NMDS-Coordinates
nmds_points <- as.data.frame(scores(NMDS, display = "sites"))
nmds_points$Sample <- rownames(nmds_points)

# Add Season
nmds_points <- nmds_points %>%
  left_join(Environmental_Data_RNA, by = c("Sample" = "SRR_ID"))

# Season as factor
nmds_points$Season <- factor(nmds_points$Season, levels = c("spring", "summer", "autumn", "winter"))

# Define Season colors
season_colors <- c("spring" = "#62b6cb", 
                   "summer" = "#ff7d00", 
                   "autumn" = "#8c1c13", 
                   "winter" = "#134074")

# Plot with Prokaryoyte Vectors
ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Season)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_polygon(data = nmds_points %>% group_by(Season) %>% slice(chull(NMDS1, NMDS2)), 
               aes(fill = Season, group = Season), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_segment(data = significant_eukaryote_vectors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  geom_text_repel(data = significant_eukaryote_vectors, aes(x = NMDS1, y = NMDS2, label = Feature),
                  size = 4, color = "black") +
  annotate("text", x = Inf, y = -Inf, label = paste("Stress =", round(NMDS$stress, 3)), 
           hjust = 1, vjust = -0.5, size = 3.5) +
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  theme_minimal() +
  labs(title = "Eukaryota", x = "NMDS1", y = "NMDS2") +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )


ggsave(
  filename = "NMDS_Euk_Tax.png",  
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
  filename = "NMDS_Euk_Tax.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)

#####################################################################################################
########################################## Prokaryota NMDS with top 15 Prokaryotes as Vector ##########################################
#####################################################################################################


####################################
### Load Required Libraries
#####################################

library(readr)
library(ggalt)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(vegan)
library(dplyr)
library(ggrepel)
library(tidyr)
library(tibble)



# Table RNA rarefied
data_rarefied_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_RNA.csv")

# Counttable (for Gernera identification)
count_data_pr2 <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_pr2.csv")
count_data_silva <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_silva.csv")

# Environmental RNA table
Environmental_Data_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metatranscriptomics.csv")
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")


#####################################
### Step 2: Extract Prokaryotes & Eukaryotes from rarefied table
#####################################

# Set SRR_ID as column
data_rarefied_RNA <- data_rarefied_RNA %>%
  rename(SRR_ID = ...1)


#####################################
### Step 3: Extract Eukaryotes Genera & map at Order
#####################################
# Nur Genera auswählen, die in data_rarefied_RNA vorkommen
eukaryote_genera_present <- unique(count_data_pr2$Genus[count_data_pr2$Genus %in% colnames(data_rarefied_RNA)])

# Filtere count_data_pr2, um nur relevante Ordnungen zu behalten
genus_to_order_pr2_filtered <- count_data_pr2 %>%
  filter(Genus %in% eukaryote_genera_present) %>%
  select(Genus, Order) %>%
  distinct()

# Daten der relevanten Genera aus der RNA-Tabelle extrahieren
data_eukaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(eukaryote_genera_present))

# Genera-Spalten zu ihrer zugehörigen Ordnung umbenennen, aber doppelte vermeiden
genus_to_order_map <- setNames(genus_to_order_pr2_filtered$Order, genus_to_order_pr2_filtered$Genus)

# ️Spalten korrekt umbenennen, indem wir eine Hilfstabelle erstellen
data_eukaryote_order <- data_eukaryote %>%
  pivot_longer(-SRR_ID, names_to = "Genus", values_to = "Abundance") %>%
  left_join(genus_to_order_pr2_filtered, by = "Genus") %>%
  group_by(SRR_ID, Order) %>%  # Gruppiere nach SRR_ID und Order (nicht Genus!)
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Order, values_from = Abundance, values_fill = 0)  # Zurück ins `wide`-Format


#####################################
### Step 3: Extract Prokaryotes Genera & map at Order
#####################################

# Nur Genera auswählen, die in data_rarefied_RNA vorkommen
prokaryote_genera_present <- unique(count_data_silva$Genus[count_data_silva$Genus %in% colnames(data_rarefied_RNA)])

# Filtere count_data_silva, um nur relevante Ordnungen zu behalten
genus_to_order_silva_filtered <- count_data_silva %>%
  filter(Genus %in% prokaryote_genera_present) %>%
  select(Genus, Order) %>%
  distinct()

# Daten der relevanten Genera aus der RNA-Tabelle extrahieren
data_prokaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(prokaryote_genera_present))

# Genera-Spalten zu ihrer zugehörigen Ordnung umbenennen, aber doppelte vermeiden
genus_to_order_map <- setNames(genus_to_order_silva_filtered$Order, genus_to_order_silva_filtered$Genus)

# Spalten korrekt umbenennen, indem wir eine Hilfstabelle erstellen
data_prokaryote_order <- data_prokaryote %>%
  pivot_longer(-SRR_ID, names_to = "Genus", values_to = "Abundance") %>%
  left_join(genus_to_order_silva_filtered, by = "Genus") %>%
  group_by(SRR_ID, Order) %>%  # Gruppiere nach SRR_ID und Order (nicht Genus!)
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Order, values_from = Abundance, values_fill = 0)  # Zurück ins `wide`-Format


#####################################
### Step 5: Remove Outlier from Metadata
#####################################

# Remove outlier samples (are already exluded from rarefied data table, son only exclude them from Env. table)
removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)

# Select only necessary Metadata columns
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season)  

# Set SRR_ID as rowname
data_prokaryote_order <- data_prokaryote_order %>% column_to_rownames(var = "SRR_ID")
data_eukaryote_order <- data_eukaryote_order %>% column_to_rownames(var = "SRR_ID")


#####################################
### Step 6: Choose optimal number of dimensions
#####################################
# Normalization of the data is neccessary!
# Although we use the rarefied table, which is already "normalized" 
# BUT! We divide pr2 and silva dataset! So new relativization is neccessary!
# Relativieren pro Probe (Zeile)
euk_data_rel <- data_eukaryote_order / rowSums(data_eukaryote_order)
# Perform NMDS with 1 to 4 dimensions and save the stress value
stress <- c()
for (i in 1:4) {
  stress <- c(stress, metaMDS(euk_data_rel, dist="bray", k = i)$stress)
}

# Plot the stress values vs dimensions   
plot(1:4, stress, xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")


#####################################
### Step 7: Calculate NMDS
#####################################
# As the final result depends on the initial random placement of the points 
# we set a seed to make the results reproducible
set.seed(183672807)

# Carry out NMDS
NMDS <- metaMDS(euk_data_rel, dist="bray", k=2, maxit=1000, trymax=100, autotransform = FALSE)
# k = number of dimensions (in this case I choose 2 Dimensions, because it refers to a stress value of a good fit ~0.15)
# trymax = maximum numbers of random starts in search of stable solution
# maxit = Maximum number of iterations in the single NMDS run
# autotransform = applies a square root transformation on the community matrix

# Plot NMDS
plot(NMDS, display = "sites")

Stress <- NMDS$stress
head(Stress)


#####################################
### Step 6: Fit Prokaryote-Taxa as Vectors
#####################################
Prok_data_rel <- data_prokaryote_order / rowSums(data_prokaryote_order)
# Extract Top 15 abundant Taxa
top_prokaryotes <- colSums(Prok_data_rel) %>%
  sort(decreasing = TRUE) %>%
  head(15)

# Filter data on top 15
prokaryote_data_top <- Prok_data_rel[, names(top_prokaryotes)]

# Fit NMDS with Prokaryotes-Taxa as Vectors
prokaryote_vectors <- envfit(NMDS, prokaryote_data_top, permutations = 999)

# Extract p-values and significant Taxa (p<0.05)
prokaryote_pvalues <- prokaryote_vectors$vectors$pvals

# Behalte nur signifikante Taxa (p < 0.05)
significant_prokaryote_vectors <- as.data.frame(scores(prokaryote_vectors, "vectors")) * ordiArrowMul(prokaryote_vectors)
significant_prokaryote_vectors$Feature <- rownames(significant_prokaryote_vectors)
significant_prokaryote_vectors <- significant_prokaryote_vectors[prokaryote_pvalues < 0.05, ]


#####################################
### Step 7: NMDS-Plot erstellen
#####################################

# Extract NMDS-Coordinates
nmds_points <- as.data.frame(scores(NMDS, display = "sites"))
nmds_points$Sample <- rownames(nmds_points)

# Add Season
nmds_points <- nmds_points %>%
  left_join(Environmental_Data_RNA, by = c("Sample" = "SRR_ID"))

# Season as factor
nmds_points$Season <- factor(nmds_points$Season, levels = c("spring", "summer", "autumn", "winter"))

# Define Season colors
season_colors <- c("spring" = "#62b6cb", 
                   "summer" = "#ff7d00", 
                   "autumn" = "#8c1c13", 
                   "winter" = "#134074")

# Plot with Prokaryoyte Vectors
ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Season)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_polygon(data = nmds_points %>% group_by(Season) %>% slice(chull(NMDS1, NMDS2)), 
               aes(fill = Season, group = Season), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_segment(data = significant_prokaryote_vectors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  geom_text_repel(data = significant_prokaryote_vectors, aes(x = NMDS1, y = NMDS2, label = Feature),
                  size = 4, color = "black") +
  annotate("text", x = Inf, y = -Inf, label = paste("Stress =", round(NMDS$stress, 3)), 
           hjust = 1, vjust = -0.5, size = 3.5) +
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  theme_minimal() +
  labs(title = "Eukaryota", x = "NMDS1", y = "NMDS2") +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )


ggsave(
  filename = "NMDS_Euk_Prok.png",  
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
  filename = "NMDS_Euk_Prok.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)


#####################################################################################################
########################################## Prokaryota NMDS with top 10 Virus as Vector ##########################################
#####################################################################################################


####################################
### Load Required Libraries
#####################################

library(readr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)



# Table RNA & DNA rarefied
data_rarefied_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_RNA.csv")
data_rarefied_refseq <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_refseq_486.csv")

# Counttable (for Gernera identification)
count_data_pr2 <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_pr2.csv")
count_data_silva <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_silva.csv")

# Environmental RNA table
Environmental_Data_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metatranscriptomics.csv")
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")


#####################################
### Step 2: Extract Prokaryotes & Eukaryotes from rarefied table
#####################################

# Set SRR_ID as column
data_rarefied_RNA <- data_rarefied_RNA %>%
  rename(SRR_ID = ...1)

data_rarefied_refseq <- data_rarefied_refseq %>% rename(SRR_ID = ...1)



#####################################
### Step 3: Extract Prokaryotes Genera & map at Order
#####################################
# Nur Genera auswählen, die in data_rarefied_RNA vorkommen
eukaryote_genera_present <- unique(count_data_pr2$Genus[count_data_pr2$Genus %in% colnames(data_rarefied_RNA)])

# Filtere count_data_pr2, um nur relevante Ordnungen zu behalten
genus_to_order_pr2_filtered <- count_data_pr2 %>%
  filter(Genus %in% eukaryote_genera_present) %>%
  select(Genus, Order) %>%
  distinct()

# Daten der relevanten Genera aus der RNA-Tabelle extrahieren
data_eukaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(eukaryote_genera_present))

# Genera-Spalten zu ihrer zugehörigen Ordnung umbenennen, aber doppelte vermeiden
genus_to_order_map <- setNames(genus_to_order_pr2_filtered$Order, genus_to_order_pr2_filtered$Genus)

# ️Spalten korrekt umbenennen, indem wir eine Hilfstabelle erstellen
data_eukaryote_order <- data_eukaryote %>%
  pivot_longer(-SRR_ID, names_to = "Genus", values_to = "Abundance") %>%
  left_join(genus_to_order_pr2_filtered, by = "Genus") %>%
  group_by(SRR_ID, Order) %>%  # Gruppiere nach SRR_ID und Order (nicht Genus!)
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Order, values_from = Abundance, values_fill = 0)  # Zurück ins `wide`-Format



#####################################
### Step 3: Remove Outlier from Metadata
#####################################

# Remove outlier samples (are already exluded from rarefied data table, son only exclude them from Env. table)
removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)

# select only neccessary Metadata columns
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season) 

# Remove outlier samples
removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
data_rarefied_refseq <- data_rarefied_refseq %>% filter(!SRR_ID %in% removed_refseq_samples)
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)



#####################################
### Step 4: Change SRR_ID Virus to Date and combine with Environmental table RNA
#####################################
# Add Date to SRR_ID of Virus to later connect witch SRR_ID RNA
data_rarefied_refseq <- left_join(data_rarefied_refseq, Environmental_Data_DNA[, c("SRR_ID", "Date")], by = "SRR_ID") %>%
  select(-SRR_ID)

# Füge SRR_ID aus Environmental_Data_RNA über Date hinzu
data_rarefied_refseq <- left_join(data_rarefied_refseq, Environmental_Data_RNA[, c("SRR_ID", "Date")], by = "Date") %>%
  select(-Date)  # Entferne Date wieder



#####################################
### Step 5: Prepare for NMDS
#####################################
# Setze SRR_ID als Zeilenname für die NMDS-Analyse
data_eukaryote_order <- data_eukaryote_order %>% column_to_rownames(var = "SRR_ID")
data_rarefied_refseq <- data_rarefied_refseq %>% column_to_rownames(var = "SRR_ID")



#####################################
### Step 6: Choose optimal number of dimensions
#####################################
# Normalization of the data is neccessary!
# Although we use the rarefied table, which is already "normalized" 
# BUT! We divide pr2 and silva dataset! So new relativization is neccessary!
# Relativieren pro Probe (Zeile)
euk_data_rel <- data_eukaryote_order / rowSums(data_eukaryote_order)
# Perform NMDS with 1 to 4 dimensions and save the stress value
stress <- c()
for (i in 1:4) {
  stress <- c(stress, metaMDS(euk_data_rel, dist="bray", k = i)$stress)
}

# Plot the stress values vs dimensions   
plot(1:4, stress, xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")



#####################################
### Step 7: Calculate NMDS
#####################################
# As the final result depends on the initial random placement of the points 
# we set a seed to make the results reproducible
set.seed(183672807)

# Carry out NMDS
NMDS <- metaMDS(euk_data_rel, dist="bray", k=2, maxit=1000, trymax=100, autotransform = FALSE)
# k = number of dimensions (in this case I choose 2 Dimensions, because it refers to a stress value of a good fit ~0.15)
# trymax = maximum numbers of random starts in search of stable solution
# maxit = Maximum number of iterations in the single NMDS run
# autotransform = applies a square root transformation on the community matrix

# Plot NMDS
plot(NMDS, display = "sites")

Stress <- NMDS$stress
head(Stress)


#####################################
### Step 4: Fit Virus-Taxa as Vectors
#####################################

# Extract Top 10 abundant Taxa
top_virus_taxa <- colSums(data_rarefied_refseq) %>%
  sort(decreasing = TRUE) %>%
  head(10)  # Wähle die 15 häufigsten Taxa

# Filter data on top 10
virus_data_top <- data_rarefied_refseq[, names(top_virus_taxa)]

# Fit NMDS mwith Virus-Taxa as Vectors
virus_vectors <- envfit(NMDS, virus_data_top, permutations = 999)

# Extract p-values and significant Taxa (p<0.05)
virus_pvalues <- virus_vectors$vectors$pvals

# Behalte nur signifikante Taxa (p < 0.05)
significant_virus_vectors <- as.data.frame(scores(virus_vectors, "vectors")) * ordiArrowMul(virus_vectors)
significant_virus_vectors$Feature <- rownames(significant_virus_vectors)
significant_virus_vectors <- significant_virus_vectors[virus_pvalues < 0.05, ]


#####################################
### Step 5: Make ggplot
#####################################
# Extract NMDS-Coordinates
nmds_points <- as.data.frame(scores(NMDS, display = "sites"))
nmds_points$Sample <- rownames(nmds_points)

# Add Season
nmds_points <- nmds_points %>%
  left_join(Environmental_Data_RNA, by = c("Sample" = "SRR_ID"))

# Season as factor
nmds_points$Season <- factor(nmds_points$Season, levels = c("spring", "summer", "autumn", "winter"))

# Define Season colors
season_colors <- c("spring" = "#62b6cb", 
                   "summer" = "#ff7d00", 
                   "autumn" = "#8c1c13", 
                   "winter" = "#134074")

# Plot with Prokaryoyte Vectors
ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Season)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_polygon(data = nmds_points %>% group_by(Season) %>% slice(chull(NMDS1, NMDS2)), 
               aes(fill = Season, group = Season), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_segment(data = significant_virus_vectors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  geom_text_repel(data = significant_virus_vectors, aes(x = NMDS1, y = NMDS2, label = Feature),
                  size = 4, color = "black") +
  annotate("text", x = Inf, y = -Inf, label = paste("Stress =", round(NMDS$stress, 3)), 
           hjust = 1, vjust = -0.5, size = 3.5) +
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  theme_minimal() +
  labs(title = "Eukaryota", x = "NMDS1", y = "NMDS2") +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )


ggsave(
  filename = "NMDS_Euk_Vir.png",  
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
  filename = "NMDS_Euk_Vir.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)

#####################################################################################################
########################################## Prokaryota NMDS with top 10 ARG as Vector ##########################################
#####################################################################################################


####################################
### Load Required Libraries
#####################################

library(readr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)



# Table RNA & DNA rarefied
data_rarefied_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_RNA.csv")

# Counttable (for Gernera identification)
count_data_silva <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_silva.csv")
count_data_pr2 <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_pr2.csv")

# Environmental RNA table
Environmental_Data_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metatranscriptomics.csv")
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")

# ARG 16S from normalized table
ARG_data <- read_csv("Master/Metatranscriptomics/Excel_lists/ARG/ARG_data_combined.csv")

#####################################
### Step 2: Extract Prokaryotes
#####################################

# Set SRR_ID as column
data_rarefied_RNA <- data_rarefied_RNA %>%
  rename(SRR_ID = ...1)



#####################################
### Step 3: Extract Eukaryotes Genera & map at Order
#####################################
# Nur Genera auswählen, die in data_rarefied_RNA vorkommen
eukaryote_genera_present <- unique(count_data_pr2$Genus[count_data_pr2$Genus %in% colnames(data_rarefied_RNA)])

# Filtere count_data_pr2, um nur relevante Ordnungen zu behalten
genus_to_order_pr2_filtered <- count_data_pr2 %>%
  filter(Genus %in% eukaryote_genera_present) %>%
  select(Genus, Order) %>%
  distinct()

# Daten der relevanten Genera aus der RNA-Tabelle extrahieren
data_eukaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(eukaryote_genera_present))

# Genera-Spalten zu ihrer zugehörigen Ordnung umbenennen, aber doppelte vermeiden
genus_to_order_map <- setNames(genus_to_order_pr2_filtered$Order, genus_to_order_pr2_filtered$Genus)

# ️Spalten korrekt umbenennen, indem wir eine Hilfstabelle erstellen
data_eukaryote_order <- data_eukaryote %>%
  pivot_longer(-SRR_ID, names_to = "Genus", values_to = "Abundance") %>%
  left_join(genus_to_order_pr2_filtered, by = "Genus") %>%
  group_by(SRR_ID, Order) %>%  # Gruppiere nach SRR_ID und Order (nicht Genus!)
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Order, values_from = Abundance, values_fill = 0)  # Zurück ins `wide`-Format


#####################################
### Step 2: Prepare ARG table
#####################################
# select only neccessary Metadata columns
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season) 

Environmental_Data_DNA <- Environmental_Data_DNA %>%
  select(SRR_ID, Date, Season)

# Konvertiere ins Long-Format
ARG_long <- ARG_data %>%
  pivot_longer(cols = -Family, names_to = "SRR_ID", values_to = "Abundance")

# Konvertiere zurück ins Wide-Format, aber mit SRR_ID als Zeilenname
ARG_transformed <- ARG_long %>%
  pivot_wider(names_from = Family, values_from = Abundance, values_fill = 0)

# Add Date to SRR_ID of Virus to later connect witch SRR_ID RNA
ARG_transformed <- left_join(ARG_transformed, Environmental_Data_DNA[, c("SRR_ID", "Date")], by = "SRR_ID") %>%
  select(-SRR_ID)

# Füge SRR_ID aus Environmental_Data_RNA über Date hinzu
ARG_transformed <- left_join(ARG_transformed, Environmental_Data_RNA[, c("SRR_ID", "Date")], by = "Date") %>%
  select(-Date)  # Entferne Date wieder

# Setze SRR_ID als Zeilennamen
ARG_matrix <- ARG_transformed %>%
  column_to_rownames("SRR_ID")


#####################################
### Step 3: Remove Outlier from Metadata
#####################################

# Remove outlier samples (are already exluded from rarefied data table, son only exclude them from Env. table)
removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)
ARG_matrix <- ARG_matrix[!rownames(ARG_matrix) %in% removed_rna_samples, ]

removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)


# select only neccessary Metadata columns
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season) 

# Remove outlier samples ARG
#removed_ARG_samples <- c("SRR1611149", "SRR1544596", "SRR1611146")
#Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_ARG_samples)

#removed_ARG_samples_RNA <- c("SRR1611150", "SRR1544599", "SRR1611147")
#Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_ARG_samples_RNA)
#data_eukaryote_order <- data_eukaryote_order %>% filter(!SRR_ID %in% removed_ARG_samples_RNA)


#####################################
### Step 4: Relativize and transform data ARG (per Sample)
#####################################
ARG_matrix <- ARG_matrix / rowSums(ARG_matrix)

# transform back to data frame
ARG_data_transformed <- as.data.frame(ARG_matrix)


#####################################
### Step 5: Prepare for NMDS
#####################################
# Setze SRR_ID als Zeilenname für die NMDS-Analyse
data_eukaryote_order <- data_eukaryote_order %>% column_to_rownames(var = "SRR_ID")


#####################################
### Step 6: Choose optimal number of dimensions
#####################################
# Normalization of the data is neccessary!
# Although we use the rarefied table, which is already "normalized" 
# BUT! We divide pr2 and silva dataset! So new relativization is neccessary!
# Relativieren pro Probe (Zeile)
euk_data_rel <- data_eukaryote_order / rowSums(data_eukaryote_order)
# Perform NMDS with 1 to 4 dimensions and save the stress value
stress <- c()
for (i in 1:4) {
  stress <- c(stress, metaMDS(euk_data_rel, dist="bray", k = i)$stress)
}

# Plot the stress values vs dimensions   
plot(1:4, stress, xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")



#####################################
### Step 7: Calculate NMDS
#####################################
# As the final result depends on the initial random placement of the points 
# we set a seed to make the results reproducible
set.seed(183672807)

# Carry out NMDS
NMDS <- metaMDS(euk_data_rel, dist="bray", k=2, maxit=1000, trymax=100, autotransform = FALSE)
# k = number of dimensions (in this case I choose 2 Dimensions, because it refers to a stress value of a good fit ~0.15)
# trymax = maximum numbers of random starts in search of stable solution
# maxit = Maximum number of iterations in the single NMDS run
# autotransform = applies a square root transformation on the community matrix

# Plot NMDS
plot(NMDS, display = "sites")

Stress <- NMDS$stress
head(Stress)


#####################################
### Step 6: ARGs as Vector
#####################################

# Calculate Top 10 abundant ARG
top_ARGs <- colSums(ARG_data_transformed) %>%
  sort(decreasing = TRUE) %>%
  head(15)

# Filtere NMDS-Data to top taxa
ARG_data_top <- ARG_data_transformed[, names(top_ARGs)]

# Fit NMDS with ARGs as Vector
ARG_vectors <- envfit(NMDS, ARG_data_top, permutations = 999)

# Extrahiere p-Werte
ARG_pvalues <- ARG_vectors$vectors$pvals

# Behalte nur signifikante ARGs (p < 0.05)
significant_ARG_vectors <- as.data.frame(scores(ARG_vectors, "vectors")) * ordiArrowMul(ARG_vectors)
significant_ARG_vectors$Feature <- rownames(significant_ARG_vectors)
significant_ARG_vectors <- significant_ARG_vectors[ARG_pvalues < 0.05, ]


#####################################
### Step 5: Make ggplot
#####################################
# Extract NMDS-Coordinates
nmds_points <- as.data.frame(scores(NMDS, display = "sites"))
nmds_points$Sample <- rownames(nmds_points)

# Add Season
nmds_points <- nmds_points %>%
  left_join(Environmental_Data_RNA, by = c("Sample" = "SRR_ID"))

# Season as factor
nmds_points$Season <- factor(nmds_points$Season, levels = c("spring", "summer", "autumn", "winter"))

# Define Season colors
season_colors <- c("spring" = "#62b6cb", 
                   "summer" = "#ff7d00", 
                   "autumn" = "#8c1c13", 
                   "winter" = "#134074")

# Plot with Prokaryoyte Vectors
ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Season)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_polygon(data = nmds_points %>% group_by(Season) %>% slice(chull(NMDS1, NMDS2)), 
               aes(fill = Season, group = Season), alpha = 0.2, color = NA) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  geom_segment(data = significant_ARG_vectors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  geom_text_repel(data = significant_ARG_vectors, aes(x = NMDS1, y = NMDS2, label = Feature),
                  size = 4, color = "black") +
  annotate("text", x = Inf, y = -Inf, label = paste("Stress =", round(NMDS$stress, 3)), 
           hjust = 1, vjust = -0.5, size = 3.5) +
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  theme_minimal() +
  labs(title = "Eukaryota", x = "NMDS1", y = "NMDS2") +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )


ggsave(
  filename = "NMDS_Euk_ARG.png",  
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
  filename = "NMDS_Euk_ARG.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)v