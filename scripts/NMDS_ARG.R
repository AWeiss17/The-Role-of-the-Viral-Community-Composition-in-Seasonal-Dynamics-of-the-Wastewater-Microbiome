########################################################################################################################
################################################# ARG NMDS  #################################
########################################################################################################################

library(readr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)




# Environmental RNA table
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")

# ARG 16S from normalized table
ARG_data <- read_csv("Master/Metatranscriptomics/Excel_lists/ARG/ARG_data_combined.csv")


# Relevant Env. variables
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten

#####################################
### Step 2: Prepare ARG table
#####################################
# Konvertiere ins Long-Format
ARG_long <- ARG_data %>%
  pivot_longer(cols = -Family, names_to = "SRR_ID", values_to = "Abundance")

# Konvertiere zurück ins Wide-Format, aber mit SRR_ID als Zeilenname
ARG_transformed <- ARG_long %>%
  pivot_wider(names_from = Family, values_from = Abundance, values_fill = 0)

# Setze SRR_ID als Zeilennamen
ARG_matrix <- ARG_transformed %>%
  column_to_rownames("SRR_ID")


#####################################
### Step 3: Remove outliers
#####################################

removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)
ARG_matrix <- ARG_matrix[!rownames(ARG_matrix) %in% removed_refseq_samples, ]

#removed_ARG_samples <- c("SRR1611149", "SRR1544596", "SRR1611146")
#Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_ARG_samples)


#####################################
### Step 4: Relativize and transform data ARG (per Sample)
#####################################
ARG_matrix <- ARG_matrix / rowSums(ARG_matrix)

# transform back to data frame
ARG_data_transformed <- as.data.frame(ARG_matrix)



#####################################
### Step 5: Calculate NMDS for Virus Data
#####################################

stress <- c()
for (i in 1:4) {
  stress <- c(stress, metaMDS(ARG_data_transformed, dist="bray", k = i)$stress)
}

# Plot the stress values vs dimensions   
plot(1:4, stress, xlab = "# of Dimensions", ylab = "Stress", main = "NMDS stress plot")

# save NMDS Stress plot
png("ARG_NMDS_StressPlot.png", width = 800, height = 600, res = 150)
plot(1:4, stress, xlab = "# of Dimensions", ylab = "Stress", main = "Virus NMDS stress plot")
dev.off()


#####################################
### Step 3: Calculate NMDS
#####################################
# As the final result depends on the initial random placement of the points 
# we set a seed to make the results reproducible
set.seed(183672807)

# Carry out NMDS
NMDS <- metaMDS(ARG_data_transformed, dist = "bray", k = 2, maxit = 1000, trymax = 100, autotransform = FALSE)

# Plot NMDS basic
plot(NMDS, display = "sites")

Stress <- NMDS$stress
head(Stress)

# NMDS-Koordinaten extrahieren
NMDS_ARG_df <- as.data.frame(scores(NMDS, display = "sites")) %>%
  rownames_to_column("SRR_ID")  # SRR_ID als Spalte hinzufügen

# Speichern der NMDS-Koordinaten als CSV
write_csv(NMDS_ARG_df, "NMDS_ARG_Coordinates.csv")
#####################################
### Step 5: Fit environmental variables onto the NMDS of refseq
#####################################

# Convert tibble to data frame if needed
Environmental_Data_DNA <- as.data.frame(Environmental_Data_DNA)

# Set `SRR_ID` as row name
rownames(Environmental_Data_DNA) <- Environmental_Data_DNA$SRR_ID

# Remove the `SRR_ID` column
Environmental_Data_DNA <- Environmental_Data_DNA[, -1]

# Make sure SRR_ID are in same order
Environmental_Data_DNA <- Environmental_Data_DNA[match(rownames(ARG_data_transformed), rownames(Environmental_Data_DNA)), ]

# Check order
identical(rownames(ARG_data_transformed), rownames(Environmental_Data_DNA)) 

# Check if there are any NA values in the dataset
anyNA(Environmental_Data_DNA)

# Fit environmental vectors and factors onto the ordination
env <- envfit(NMDS, Environmental_Data_DNA, permutations = 999)

# Plot environmental vectors
plot(NMDS, display = "sites")
plot(env, p.max = 0.05)



#####################################
### Step 6: Extract only significant environmental facors
#####################################

# Extract p-values from environmental factors
env_pvalues <- env$vectors$pvals
print(env$vectors$pvals)

# Filtere nur signifikante Vektoren (p < 0.05)
significant_env_vectors <- as.data.frame(scores(env, "vectors")) * ordiArrowMul(env)
significant_env_vectors$Feature <- rownames(significant_env_vectors)
significant_env_vectors <- significant_env_vectors[env_pvalues < 0.05, ]

# Erstelle eine neue Spalte mit den gewünschten Bezeichnungen
label_mapping <- c("pH_online" = "pH","PO4_online_mg.L" = "Phosphate", "NO3_online_mg.L" = "Nitrate", "Temp_manual" = "Water temperature", "Temp_Air" = "Air temperature")

# Neue Labels in die `significant_env_vectors`-Tabelle übernehmen
significant_env_vectors$Feature_renamed <- label_mapping[significant_env_vectors$Feature]

# Falls ein Feature nicht in `label_mapping` ist, behalte den Originalnamen
significant_env_vectors$Feature_renamed[is.na(significant_env_vectors$Feature_renamed)] <- significant_env_vectors$Feature


#####################################
### Step 7: Make ggplot
#####################################
# Extact NMDS-Coordinates
nmds_points <- as.data.frame(scores(NMDS, display = "sites"))
nmds_points$Sample <- rownames(nmds_points)

# Add Season
nmds_points$Season <- Environmental_Data_DNA$Season[match(rownames(nmds_points), rownames(Environmental_Data_DNA))]

# Season as factor
nmds_points$Season <- as.factor(nmds_points$Season)

# set order within legend
nmds_points$Season <- factor(nmds_points$Season, levels = c("spring", "summer", "autumn", "winter"))

# Define Season colors
season_colors <- c("spring" = "#62b6cb", 
                   "summer" = "#ff7d00", 
                   "autumn" = "#8c1c13", 
                   "winter" = "#134074")


# Extract environmental factors from envit
env_vectors <- as.data.frame(scores(env, "vectors")) * ordiArrowMul(env)
env_vectors$Feature <- rownames(env_vectors)

# Plot with environmental factors
ggplot(nmds_points, aes(x = NMDS1, y = NMDS2, color = Season)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_polygon(data = nmds_points %>% group_by(Season) %>% slice(chull(NMDS1, NMDS2)), 
               aes(fill = Season, group = Season), alpha = 0.2, color = NA) +
  # Dezente graue Linien für den Nullpunkt
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") +
  geom_vline(xintercept = 0, color = "gray", linetype = "dashed") +
  
  # Nur signifikante Vektoren für Umweltfaktoren
  geom_segment(data = significant_env_vectors, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), color = "black") +
  
  # Beschriftung der signifikanten Umweltfaktoren
  geom_text_repel(data = significant_env_vectors, aes(x = NMDS1, y = NMDS2, label = Feature_renamed),
                  size = 4, color = "black") +
  
  # Stress-Wert EXAKT unten rechts setzen (ohne den Plot zu verkleinern!)
  annotate("text", 
           x = Inf,  # Ganz rechts
           y = -Inf,  # Ganz unten
           label = paste("Stress value:", round(NMDS$stress, 3)), 
           hjust = 1, vjust = -0.5, size = 3.5) +  # Falls zu nah am Rand, vjust anpassen
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  theme_minimal() +
  labs(title = "ARG", x = "NMDS1", y = "NMDS2") +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )

ggsave(
  filename = "NMDS_ARG_Env.png",  
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
  filename = "NMDS_ARG_Env.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)


########################################################################################################################
################################################# ARG NMDS with Top 10 ARG Families  #################################
########################################################################################################################

library(readr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)




# Environmental RNA table
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")

# ARG 16S from normalized table
ARG_data <- read_csv("Master/Metatranscriptomics/Excel_lists/ARG/ARG_data_combined.csv")


# Relevant Env. variables
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten

#####################################
### Step 2: Prepare ARG table
#####################################

# Konvertiere ins Long-Format
ARG_long <- ARG_data %>%
  pivot_longer(cols = -Family, names_to = "SRR_ID", values_to = "Abundance")

# Konvertiere zurück ins Wide-Format, aber mit SRR_ID als Zeilenname
ARG_transformed <- ARG_long %>%
  pivot_wider(names_from = Family, values_from = Abundance, values_fill = 0)

# Setze SRR_ID als Zeilennamen
ARG_matrix <- ARG_transformed %>%
  column_to_rownames("SRR_ID")


#####################################
### Step 3: Remove outliers
#####################################

removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)
ARG_matrix <- ARG_matrix[!rownames(ARG_matrix) %in% removed_refseq_samples, ]

#removed_ARG_samples <- c("SRR1611149", "SRR1544596", "SRR1611146")
#Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_ARG_samples)


#####################################
### Step 4: Relativize and transform data ARG (per Sample)
#####################################
ARG_matrix <- ARG_matrix / rowSums(ARG_matrix)

# transform back to data frame
ARG_data_transformed <- as.data.frame(ARG_matrix)



#####################################
### Step 5: Calculate NMDS for Virus Data
#####################################

stress <- c()
for (i in 1:4) {
  stress <- c(stress, metaMDS(ARG_data_transformed, dist="bray", k = i)$stress)
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
NMDS <- metaMDS(ARG_data_transformed, dist = "bray", k = 2, maxit = 1000, trymax = 100, autotransform = FALSE)

# Plot NMDS basic
plot(NMDS, display = "sites")

Stress <- NMDS$stress
head(Stress)


#####################################
### Step 4: Fit Virus-Taxa as Vectors
#####################################

# Extract Top 10 abundant Taxa
ARG_data_transformed_top <- colSums(ARG_data_transformed) %>%
  sort(decreasing = TRUE) %>%
  head(10)  # Wähle die 15 häufigsten Taxa

# Filter data on top 10
ARG_data_top <- ARG_data_transformed[, names(ARG_data_transformed_top)]

# Fit NMDS mwith Virus-Taxa as Vectors
ARG_vectors <- envfit(NMDS, ARG_data_top, permutations = 999)

# Extract p-values and significant Taxa (p<0.05)
ARG_pvalues <- ARG_vectors$vectors$pvals

# Behalte nur signifikante Taxa (p < 0.05)
significant_ARG_vectors <- as.data.frame(scores(ARG_vectors, "vectors")) * ordiArrowMul(ARG_vectors)
significant_ARG_vectors$Feature <- rownames(significant_ARG_vectors)
significant_ARG_vectors <- significant_ARG_vectors[ARG_pvalues < 0.05, ]



#####################################
### Step 5: Make ggplot
#####################################
# Extact NMDS-Coordinates
nmds_points <- as.data.frame(scores(NMDS, display = "sites"))
nmds_points$Sample <- rownames(nmds_points)

# Add Season
nmds_points <- nmds_points %>%
  left_join(Environmental_Data_DNA, by = c("Sample" = "SRR_ID"))

# Season as factor
nmds_points$Season <- as.factor(nmds_points$Season)

# set order within legend
nmds_points$Season <- factor(nmds_points$Season, levels = c("spring", "summer", "autumn", "winter"))

# Define Season colors
season_colors <- c("spring" = "#62b6cb", 
                   "summer" = "#ff7d00", 
                   "autumn" = "#8c1c13", 
                   "winter" = "#134074")


# Plot with viral taxa
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
  annotate("text", 
           x = Inf,  
           y = -Inf,  
           label = paste("Stress value:", round(NMDS$stress, 3)), 
           hjust = 1, vjust = -0.5, size = 3.5) +  
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  theme_minimal() +
  labs(title = "ARG", x = "NMDS1", y = "NMDS2") +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )

ggsave(
  filename = "NMDS_ARG_Tax.png",  
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
  filename = "NMDS_ARG_Taxa.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)


#########################################################################################################################
################################################# Viral NMDS with top 15 prokaryotes as Vector #################################
#########################################################################################################################

library(readr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)



# Table RNA & DNA rarefied
data_rarefied_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_RNA.csv")
# ARG 16S from normalized table
ARG_data <- read_csv("Master/Metatranscriptomics/Excel_lists/ARG/ARG_data_combined.csv")

# Counttable (for Gernera identification)
count_data_silva <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_silva.csv")

# Environmental RNA table
Environmental_Data_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metatranscriptomics.csv")
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")


#####################################
### Step 2: Extract Prokaryotes & Eukaryotes from rarefied table
#####################################

# Set SRR_ID as column
data_rarefied_RNA <- data_rarefied_RNA %>% rename(SRR_ID = ...1)


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

# Set SRR_ID as rowname and convert to Matrix
Prok_data_matrix <- data_prokaryote_order %>%
  column_to_rownames(var = "SRR_ID") %>%
  as.matrix()

# Normalization of the data is neccessary!
# Although we use the rarefied table, which is already "normalized" 
# BUT! We divide pr2 and silva dataset! So new relativization is neccessary!
# Relativieren pro Probe (Zeile)
Prok_data_rel <- Prok_data_matrix / rowSums(Prok_data_matrix)

# Matrix zurück in DataFrame
Prok_data_table <- as.data.frame(Prok_data_rel)

# SRR_ID aus rownames wieder als eigene Spalte
Prok_data_table <- Prok_data_table %>%
  rownames_to_column(var = "SRR_ID")


#####################################
### Step 2: Prepare ARG table
#####################################

# Konvertiere ins Long-Format
ARG_long <- ARG_data %>%
  pivot_longer(cols = -Family, names_to = "SRR_ID", values_to = "Abundance")

# Konvertiere zurück ins Wide-Format, aber mit SRR_ID als Zeilenname
ARG_transformed <- ARG_long %>%
  pivot_wider(names_from = Family, values_from = Abundance, values_fill = 0)

# Setze SRR_ID als Zeilennamen
ARG_matrix <- ARG_transformed %>%
  column_to_rownames("SRR_ID")


#####################################
### Step 3: Remove outliers
#####################################
# Remove outlier samples (are already exluded from rarefied data table, son only exclude them from Env. table)
removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)

removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)
ARG_matrix <- ARG_matrix[!rownames(ARG_matrix) %in% removed_refseq_samples, ]

# select only neccessary Metadata columns
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season) 

# Remove outlier samples ARG
#removed_ARG_samples <- c("SRR1611149", "SRR1544596", "SRR1611146")
#Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_ARG_samples)

#removed_ARG_samples_RNA <- c("SRR1611150", "SRR1544599", "SRR1611147")
#Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_ARG_samples_RNA)
#Prok_data_table <- Prok_data_table %>% filter(!SRR_ID %in% removed_ARG_samples_RNA)




#####################################
### Step 4: Relativize and transform data ARG (per Sample)
#####################################
ARG_matrix <- ARG_matrix / rowSums(ARG_matrix)

# transform back to data frame
ARG_data_transformed <- as.data.frame(ARG_matrix)



#####################################
### Step 5: Calculate NMDS for Virus Data
#####################################

stress <- c()
for (i in 1:4) {
  stress <- c(stress, metaMDS(ARG_data_transformed, dist="bray", k = i)$stress)
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
NMDS <- metaMDS(ARG_data_transformed, dist = "bray", k = 2, maxit = 1000, trymax = 100, autotransform = FALSE)

# Plot NMDS basic
plot(NMDS, display = "sites")

Stress <- NMDS$stress
head(Stress)


#####################################
### Step 4: Change SRR_ID Prokaryotes to Date and combine with Environmental table DNA
#####################################
# Add Date to SRR_ID of Virus to later connect witch SRR_ID RNA
Prok_data_table <- left_join(Prok_data_table, Environmental_Data_RNA[, c("SRR_ID", "Date")], by = "SRR_ID") %>%
  select(-SRR_ID)

# Füge SRR_ID aus Environmental_Data_DNA über Date hinzu
Prok_data_table <- left_join(Prok_data_table, Environmental_Data_DNA[, c("SRR_ID", "Date")], by = "Date") %>%
  select(-Date)  # Entferne Date wieder

# Setze SRR_ID als Zeilennamen
Prok_data_table <- Prok_data_table %>%
  column_to_rownames("SRR_ID")



#####################################
### Step 6: Prokaryoten als Vektoren anpassen
#####################################
# Wähle die 15 häufigsten Prokaryoten
top_prokaryotes <- colSums(Prok_data_table) %>%
  sort(decreasing = TRUE) %>%
  head(15)

# Filtere NMDS-Daten auf Top 15 Prokaryoten
prokaryote_data_top <- Prok_data_table[, names(top_prokaryotes)]

# Fit NMDS mit Prokaryoten als Vektoren
prokaryote_vectors <- envfit(NMDS, prokaryote_data_top, permutations = 999)

# Extrahiere p-Werte
prokaryote_pvalues <- prokaryote_vectors$vectors$pvals

# Behalte nur signifikante Taxa (p < 0.05)
significant_prokaryote_vectors <- as.data.frame(scores(prokaryote_vectors, "vectors")) * ordiArrowMul(prokaryote_vectors)
significant_prokaryote_vectors$Feature <- rownames(significant_prokaryote_vectors)
significant_prokaryote_vectors <- significant_prokaryote_vectors[prokaryote_pvalues < 0.05, ]



#####################################
### Step 7: NMDS-Plot erstellen
#####################################
# NMDS-Koordinaten extrahieren
nmds_points <- as.data.frame(scores(NMDS, display = "sites"))
nmds_points$Sample <- rownames(nmds_points)

# Seasons hinzufügen
nmds_points <- nmds_points %>%
  left_join(Environmental_Data_DNA, by = c("Sample" = "SRR_ID"))

# Sicherstellen, dass Season ein Faktor ist
nmds_points$Season <- factor(nmds_points$Season, levels = c("spring", "summer", "autumn", "winter"))

# Farben für Seasons
season_colors <- c("spring" = "#62b6cb",
                   "summer" = "#ff7d00",
                   "autumn" = "#8c1c13", 
                   "winter" = "#134074")

# Plot
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
  annotate("text", x = Inf, y = -Inf, label = paste("Stress value:", round(NMDS$stress, 3)), 
           hjust = 1, vjust = -0.5, size = 3.5) +
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  theme_minimal() +
  labs(title = "ARG", x = "NMDS1", y = "NMDS2") +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )


ggsave(
  filename = "NMDS_ARG_Prok.png",  
  plot = last_plot(),  
  device = "png",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)



#########################################################################################################################
################################################# Viral NMDS with top 15 eukaryotes as Vector #################################
#########################################################################################################################

library(readr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)



# Table RNA & DNA rarefied
data_rarefied_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_RNA.csv")

# ARG 16S from normalized table
ARG_data <- read_csv("Master/Metatranscriptomics/Excel_lists/ARG/ARG_data_combined.csv")

# Counttable (for Gernera identification)
count_data_pr2 <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/counttable_metatranscr01_genus_new_renamed_pr2.csv")

# Environmental RNA table
Environmental_Data_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metatranscriptomics.csv")
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")


#####################################
### Step 2: Extract Eukaryotes from rarefied table
#####################################

# Set SRR_ID as column
data_rarefied_RNA <- data_rarefied_RNA %>% rename(SRR_ID = ...1)


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

# Set SRR_ID as rowname and convert to Matrix
Euk_data_matrix <- data_eukaryote_order %>%
  column_to_rownames(var = "SRR_ID") %>%
  as.matrix()

# Normalization of the data is neccessary!
# Although we use the rarefied table, which is already "normalized" 
# BUT! We divide pr2 and silva dataset! So new relativization is neccessary!
# Relativieren pro Probe (Zeile)
Euk_data_rel <- Euk_data_matrix / rowSums(Euk_data_matrix)

# Matrix zurück in DataFrame
Euk_data_table <- as.data.frame(Euk_data_rel)

# SRR_ID aus rownames wieder als eigene Spalte
Euk_data_table <- Euk_data_table %>%
  rownames_to_column(var = "SRR_ID")


#####################################
### Step 2: Prepare ARG table
#####################################

# Konvertiere ins Long-Format
ARG_long <- ARG_data %>%
  pivot_longer(cols = -Family, names_to = "SRR_ID", values_to = "Abundance")

# Konvertiere zurück ins Wide-Format, aber mit SRR_ID als Zeilenname
ARG_transformed <- ARG_long %>%
  pivot_wider(names_from = Family, values_from = Abundance, values_fill = 0)

# Setze SRR_ID als Zeilennamen
ARG_matrix <- ARG_transformed %>%
  column_to_rownames("SRR_ID")


#####################################
### Step 3: Remove outliers
#####################################
# Remove outlier samples (are already exluded from rarefied data table, son only exclude them from Env. table)
removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)

removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)
ARG_matrix <- ARG_matrix[!rownames(ARG_matrix) %in% removed_refseq_samples, ]

# select only neccessary Metadata columns
Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season) 

# Remove outlier samples ARG
#removed_ARG_samples <- c("SRR1611149", "SRR1544596", "SRR1611146")
#Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_ARG_samples)

#removed_ARG_samples_RNA <- c("SRR1611150", "SRR1544599", "SRR1611147")
#Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_ARG_samples_RNA)
#Euk_data_table <- Euk_data_table %>% filter(!SRR_ID %in% removed_ARG_samples_RNA)



#####################################
### Step 4: Relativize and transform data ARG (per Sample)
#####################################
ARG_matrix <- ARG_matrix / rowSums(ARG_matrix)

# transform back to data frame
ARG_data_transformed <- as.data.frame(ARG_matrix)



#####################################
### Step 5: Calculate NMDS for Virus Data
#####################################

stress <- c()
for (i in 1:4) {
  stress <- c(stress, metaMDS(ARG_data_transformed, dist="bray", k = i)$stress)
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
NMDS <- metaMDS(ARG_data_transformed, dist = "bray", k = 2, maxit = 1000, trymax = 100, autotransform = FALSE)

# Plot NMDS basic
plot(NMDS, display = "sites")

Stress <- NMDS$stress
head(Stress)


#####################################
### Step 4: Change SRR_ID Prokaryotes to Date and combine with Environmental table DNA
#####################################
# Add Date to SRR_ID of Virus to later connect witch SRR_ID RNA
Euk_data_table <- left_join(Euk_data_table, Environmental_Data_RNA[, c("SRR_ID", "Date")], by = "SRR_ID") %>%
  select(-SRR_ID)

# Füge SRR_ID aus Environmental_Data_DNA über Date hinzu
Euk_data_table <- left_join(Euk_data_table, Environmental_Data_DNA[, c("SRR_ID", "Date")], by = "Date") %>%
  select(-Date)  # Entferne Date wieder

# Setze SRR_ID als Zeilennamen
Euk_data_table <- Euk_data_table %>%
  column_to_rownames("SRR_ID")



#####################################
### Step 6: Prokaryoten als Vektoren anpassen
#####################################
# Extract Top 15 abundant Taxa
top_eukaryotes <- colSums(Euk_data_table) %>%
  sort(decreasing = TRUE) %>%
  head(15)

# Filter data on top 15
eukaryote_data_top <- Euk_data_table[, names(top_eukaryotes)]

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
# NMDS-Koordinaten extrahieren
nmds_points <- as.data.frame(scores(NMDS, display = "sites"))
nmds_points$Sample <- rownames(nmds_points)

# Seasons hinzufügen
nmds_points <- nmds_points %>%
  left_join(Environmental_Data_DNA, by = c("Sample" = "SRR_ID"))

# Sicherstellen, dass Season ein Faktor ist
nmds_points$Season <- factor(nmds_points$Season, levels = c("spring", "summer", "autumn", "winter"))

# Farben für Seasons
season_colors <- c("spring" = "#62b6cb",
                   "summer" = "#ff7d00",
                   "autumn" = "#8c1c13", 
                   "winter" = "#134074")

# Plot
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
  annotate("text", x = Inf, y = -Inf, label = paste("Stress value:", round(NMDS$stress, 3)), 
           hjust = 1, vjust = -0.5, size = 3.5) +
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  theme_minimal() +
  labs(title = "ARG", x = "NMDS1", y = "NMDS2") +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )


ggsave(
  filename = "NMDS_ARG_Euk.png",  
  plot = last_plot(),  
  device = "png",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)



########################################################################################################################
################################################# ARG NMDS with Top 10 Virus Families  #################################
########################################################################################################################

library(readr)
library(ggplot2)
library(ggrepel)
library(vegan)
library(dplyr)
library(tidyr)
library(tibble)




# Environmental RNA table
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")

# ARG 16S from normalized table
ARG_data <- read_csv("Master/Metatranscriptomics/Excel_lists/ARG/ARG_data_combined.csv")
data_rarefied_refseq <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_refseq_486.csv")

# Relevant Env. variables
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  select(SRR_ID, Date, Season)  


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

# Setze SRR_ID als Zeilennamen
ARG_matrix <- ARG_transformed %>%
  column_to_rownames("SRR_ID")


#####################################
### Step 3: Remove outliers
#####################################

removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)
ARG_matrix <- ARG_matrix[!rownames(ARG_matrix) %in% removed_refseq_samples, ]
data_rarefied_refseq <- data_rarefied_refseq %>% filter(!SRR_ID %in% removed_refseq_samples)

#removed_ARG_samples <- c("SRR1611149", "SRR1544596", "SRR1611146")
#Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_ARG_samples)
#data_rarefied_refseq <- data_rarefied_refseq %>% filter(!SRR_ID %in% removed_ARG_samples)


# Setze SRR_ID als Zeilenname für die NMDS-Analyse
data_rarefied_refseq <- data_rarefied_refseq %>% column_to_rownames(var = "SRR_ID")

#####################################
### Step 4: Relativize and transform data ARG (per Sample)
#####################################
ARG_matrix <- ARG_matrix / rowSums(ARG_matrix)

# transform back to data frame
ARG_data_transformed <- as.data.frame(ARG_matrix)



#####################################
### Step 5: Calculate NMDS for Virus Data
#####################################

stress <- c()
for (i in 1:4) {
  stress <- c(stress, metaMDS(ARG_data_transformed, dist="bray", k = i)$stress)
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
NMDS <- metaMDS(ARG_data_transformed, dist = "bray", k = 2, maxit = 1000, trymax = 100, autotransform = FALSE)

# Plot NMDS basic
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

print(virus_vectors)
# Extract p-values and significant Taxa (p<0.05)
virus_pvalues <- virus_vectors$vectors$pvals

# Behalte nur signifikante Taxa (p < 0.05)
significant_virus_vectors <- as.data.frame(scores(virus_vectors, "vectors")) * ordiArrowMul(virus_vectors)
significant_virus_vectors$Feature <- rownames(significant_virus_vectors)
significant_virus_vectors <- significant_virus_vectors[virus_pvalues < 0.05, ]



#####################################
### Step 5: Make ggplot
#####################################
# Extact NMDS-Coordinates
nmds_points <- as.data.frame(scores(NMDS, display = "sites"))
nmds_points$Sample <- rownames(nmds_points)

# Add Season
nmds_points <- nmds_points %>%
  left_join(Environmental_Data_DNA, by = c("Sample" = "SRR_ID"))

# Season as factor
nmds_points$Season <- as.factor(nmds_points$Season)

# set order within legend
nmds_points$Season <- factor(nmds_points$Season, levels = c("spring", "summer", "autumn", "winter"))

# Define Season colors
season_colors <- c("spring" = "#62b6cb", 
                   "summer" = "#ff7d00", 
                   "autumn" = "#8c1c13", 
                   "winter" = "#134074")


# Plot with viral taxa
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
  annotate("text", 
           x = Inf,  
           y = -Inf,  
           label = paste("Stress value:", round(NMDS$stress, 3)), 
           hjust = 1, vjust = -0.5, size = 3.5) +  
  scale_color_manual(values = season_colors) +
  scale_fill_manual(values = season_colors) +
  theme_minimal() +
  labs(title = "ARG", x = "NMDS1", y = "NMDS2") +
  theme(
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(hjust = 0, size = 14, face = "bold")
  )

ggsave(
  filename = "NMDS_ARG_Vir.png",  
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
  filename = "NMDS_Virus_Taxa.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)
