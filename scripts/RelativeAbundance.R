#####################################################################################################
########################################## Viruses ##################################################
#####################################################################################################


####################################
### Load Required Libraries
#####################################
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(car)



#####################################
### Step 1: Load rarefied counttables & environmental tables
#####################################

data_rarefied_refseq <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_refseq_486.csv")

Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")


# Set SRR_ID as column
data_rarefied_refseq <- data_rarefied_refseq %>%
  rename(SRR_ID = ...1)


# Remove outlier samples
removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
data_rarefied_refseq <- data_rarefied_refseq %>% filter(!SRR_ID %in% removed_refseq_samples)
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)


# Summe aller Reads über alle Proben und alle Taxa:
total_reads_refseq <- sum(data_rarefied_refseq[,-1], na.rm = TRUE)

print(total_reads_refseq)

#####################################
### Step 2: Set Row Names Refseq & environmental table
#####################################

data_rarefied_refseq <- as.data.frame(data_rarefied_refseq)
rownames(data_rarefied_refseq) <- data_rarefied_refseq$SRR_ID
data_rarefied_refseq <- data_rarefied_refseq[, -1]


#####################################
### Step 3: Relativize data
#####################################
# Compute relative abundance only for taxonomic counts
data_relative_refseq <- data_rarefied_refseq / rowSums(data_rarefied_refseq)



#####################################
### Step 4: Convert to long format
#####################################
# Convert data to long format for ggplot2
data_long <- data_relative_refseq %>%
  tibble::rownames_to_column(var = "SRR_ID") %>%
  pivot_longer(-SRR_ID, names_to = "Taxon", values_to = "Abundance")  

# Summarize top 15 most abundant taxa across all samples
top_taxa <- colSums(data_relative_refseq) %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  head(10)

print(top_taxa)

# Group non-top taxa into "Other"
data_long <- data_long %>%
  mutate(Taxon = ifelse(Taxon %in% top_taxa, Taxon, "Other"))

# Merge with environmental metadata 
data_long <- left_join(data_long, Environmental_Data_DNA %>% select(SRR_ID, Date, Season), by = "SRR_ID")

# Convert Date into correct format
data_long$Date <- as.Date(data_long$Date, format = "%Y-%m-%d")




#####################################
### Step 5: Plot relative abundances over Seasons
#####################################
data_long$Taxon <- factor(data_long$Taxon, levels = c("Other", setdiff(unique(data_long$Taxon), "Other")))

# Set order of Season
data_long$Season <- factor(data_long$Season, levels = c("winter", "spring", "summer", "autumn"))


custom_colors <- c(
  "Ackermannviridae" = "#113e90",  
  "Autographiviridae" = "#417de9",  
  #"Casjensviridae" = "#ffed6f",  
  "Demerecviridae" = "#bcd1f7",  
  #"Herelleviridae" = "#cce3de,  
  "Inoviridae" = "#577590", 
  "Iridoviridae" = "#ffcdb2",  
  "Kyanoviridae" = "#fcf6bd",  
  #"Marseilleviridae" = "#ffb4a2",  
  #"Mesyanzhinovviridae" = "#e5989b", 
  "Mimiviridae" = "#e56b6f", 
  #"Pachyviridae" = "#6d6875", 
  "Peduoviridae" = "#f9c74f", 
  "Schitoviridae" = "#b56576", 
  "Straboviridae" = "#43aa8b", 
  "Other" = "#403d39"    
)


#####################################
### Step 6: Plot relative abundances over all Dates
#####################################
# Coverte `Date` into correct format and sort Date
data_long$Date <- as.Date(data_long$Date, format = "%Y-%m-%d")
data_long <- data_long %>% arrange(Date)  

# Find first date in month
first_dates <- data_long %>%
  group_by(format(Date, "%Y-%m")) %>%
  summarise(first_date = min(Date)) %>%
  pull(first_date)

# create plot
ggplot(data_long, aes(x = as.factor(Date), y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) + 
  scale_x_discrete(labels = format(as.Date(unique(data_long$Date)), "%Y-%m-%d")) + 
  scale_fill_manual(values = custom_colors, name = "Taxon") +
  labs(
    title = "Virus",
    x = "Date",
    y = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave(
  filename = "Rel_Abundance_Virus.png",  
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
  filename = "Rel_Abundance_Virus.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)

#####################################
###### Statistical testing
#####################################
# Normal distribution?
shapiro.test(data_long$Abundance) 

#### not different between seasons
## Peduoviridae
Peduo_data <- data_long %>%
  filter(Taxon == "Peduoviridae")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Peduo_data)

## Inoviridae
Ino_data <- data_long %>%
  filter(Taxon == "Inoviridae")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Ino_data)


## Inoviridae
Irido_data <- data_long %>%
  filter(Taxon == "Iridoviridae")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Irido_data)


#### different between seasons
## Kyanoviridae
kyano_data <- data_long %>%
  filter(Taxon == "Kyanoviridae")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = kyano_data)

# Post-hoc: Pairwise Wilcoxon
pairwise.wilcox.test(kyano_data$Abundance, kyano_data$Season, p.adjust.method = "BH")
aggregate(Abundance ~ Season, data = kyano_data, FUN = mean)

## Autographiviridae
Auto_data <- data_long %>%
  filter(Taxon == "Autographiviridae")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Auto_data)
aggregate(Abundance ~ Season, data = Auto_data, FUN = mean)

## Mimiviridae
Mimi_data <- data_long %>%
  filter(Taxon == "Mimiviridae")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Mimi_data)

## Ackermannviridae
Acker_data <- data_long %>%
  filter(Taxon == "Ackermannviridae")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Acker_data)



#####################################################################################################
########################################## Eukaryotess ##############################################
#####################################################################################################


####################################
### Load Required Libraries
#####################################
library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)


#####################################
### Step 1: Load rarefied counttables & environmental tables
#####################################

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
# Select only those genera, which occur in data_rarefied_RNA 
prokaryote_genera_present <- unique(count_data_silva$Genus[count_data_silva$Genus %in% colnames(data_rarefied_RNA)])

# Filter count_data_silva, for relevant orders
genus_to_order_silva_filtered <- count_data_silva %>%
  filter(Genus %in% prokaryote_genera_present) %>%
  select(Genus, Order) %>%
  distinct()

# extract relevant orders from RNA table
data_prokaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(prokaryote_genera_present))

# rename Genera-columns to specific order
genus_to_order_map <- setNames(genus_to_order_silva_filtered$Order, genus_to_order_silva_filtered$Genus)

# rename columns
data_prokaryote_order <- data_prokaryote %>%
  pivot_longer(-SRR_ID, names_to = "Genus", values_to = "Abundance") %>%
  left_join(genus_to_order_silva_filtered, by = "Genus") %>%
  group_by(SRR_ID, Order) %>%  # Group after SRR_ID and Order (not Genus!)
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Order, values_from = Abundance, values_fill = 0)  


#####################################
### Step 4: Extract Eukaryotes Genera & map at Order
#####################################
# select only genera, which occur in data_rarefied_RNA vorkommen
eukaryote_genera_present <- unique(count_data_pr2$Genus[count_data_pr2$Genus %in% colnames(data_rarefied_RNA)])

# Filter count_data_pr2, for relevant order
genus_to_order_pr2_filtered <- count_data_pr2 %>%
  filter(Genus %in% eukaryote_genera_present) %>%
  select(Genus, Order) %>%
  distinct()

# extract relevant genera from RNA table
data_eukaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(eukaryote_genera_present))

# rename genera columns to specific order
genus_to_order_map <- setNames(genus_to_order_pr2_filtered$Order, genus_to_order_pr2_filtered$Genus)

# ️rename column
data_eukaryote_order <- data_eukaryote %>%
  pivot_longer(-SRR_ID, names_to = "Genus", values_to = "Abundance") %>%
  left_join(genus_to_order_pr2_filtered, by = "Genus") %>%
  group_by(SRR_ID, Order) %>%  
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Order, values_from = Abundance, values_fill = 0)  




#####################################
### Step 5: Remove Outlier from Metadata
#####################################

# Remove outlier samples (are already exluded from rarefied data table, so only exclude them from Env. table)
removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)

removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
data_rarefied_refseq <- data_rarefied_refseq %>% filter(!SRR_ID %in% removed_refseq_samples)
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)





#####################################
### Step 6: Set Row Names Refseq & environmental table
#####################################

data_eukaryote_order <- as.data.frame(data_eukaryote_order)
rownames(data_eukaryote_order) <- data_eukaryote_order$SRR_ID
data_eukaryote_order <- data_eukaryote_order[, -1]


#####################################
### Step 7: Relativize data
#####################################
# Compute relative abundance only for taxonomic counts
data_eukaryote_order <- data_eukaryote_order / rowSums(data_eukaryote_order)




#####################################
### Step 8: Convert to long format
#####################################
# Convert data to long format for ggplot2
data_long <- data_eukaryote_order %>%
  tibble::rownames_to_column(var = "SRR_ID") %>%
  pivot_longer(-SRR_ID, names_to = "Taxon", values_to = "Abundance")  

# Summarize top 15 most abundant taxa across all samples
top_taxa <- colSums(data_eukaryote_order) %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  head(15)

print(top_taxa)

# Group non-top taxa into "Other"
data_long <- data_long %>%
  mutate(Taxon = ifelse(Taxon %in% top_taxa, Taxon, "Other"))

# Merge with environmental metadata 
data_long <- left_join(data_long, Environmental_Data_RNA %>% select(SRR_ID, Date, Season), by = "SRR_ID")


# Konvertiere Date ins richtige Format
data_long$Date <- as.Date(data_long$Date, format = "%Y-%m-%d")




#####################################
### Step 9: Plot relative abundances over Seasons
#####################################
data_long$Taxon <- factor(data_long$Taxon, levels = c("Other", setdiff(unique(data_long$Taxon), "Other")))

# Set order of seasons
data_long$Season <- factor(data_long$Season, levels = c("winter", "spring", "summer", "autumn"))


custom_colors <- c(
  "Vannellida" = "#ff8fa3",  
  "Peritrichia" = "#ff4d6d",  
  "Saccharomycotina" = "#a4133c",  
  "Kinetoplastida" = "#590d22",  
  "Pezizomycotina" = "#80ffdb",  
  "Himatismenida" = "#e7bc91",  
  "Dactylopodida" = "#5390d9",  
  "Microsporidiomycotina" = "#6930c3",  
  "Trebouxiophyceae_X" = "#3c096c",  
  "Rotifera_X" = "#9d4edd", 
  "Prostomatea_X" = "#e7c6ff", 
  "Fungi_XX" = "#606c38", 
  "Plagiopylea_X" = "#ff8500", 
  "Cryomonadida" = "#52b788", 
  "Armophorea_X" = "#ffd000",
  "Nolandida" = "#ff4800",
  "Arcellinida_X" = "#ff8500",
  "Dermamoebida" = "#ffaa00",
  "Chlamydomonadales" = "#ffd000",
  "Chytridiomycotina" = "#ffea00",
  "Kickxellomycotina" = "#e7bc91",
  "Ciliophora_XX" = "#bc8a5f",
  "Syndiniales_X" = "#6f4518",
  "Euplotia" = "#748cab",
  "Haptoria" = "#3e5c76",
  "Chrysophyceae_X" = "#1d2d44",
  "Breviatidea" = "#a5be00",
  "Suctoria" = "#9d8189",
  "Granofilosea_X" = "#706677",
  "Euglenida" = "#a53860",
  "Other" = "#403d39"    
)



#####################################
### Step 10: Plot relative abundances over all Dates
#####################################
# convert `Date` into correct format and sort after date
data_long$Date <- as.Date(data_long$Date, format = "%Y-%m-%d")
data_long <- data_long %>% arrange(Date)  

# Find first date in month
first_dates <- data_long %>%
  group_by(format(Date, "%Y-%m-%d")) %>%
  summarise(first_date = min(Date)) %>%
  pull(first_date)


ggplot(data_long, aes(x = as.factor(Date), y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +  
  scale_x_discrete(labels = format(as.Date(unique(data_long$Date)), "%Y-%m-%d")) +  
  scale_fill_manual(values = custom_colors, name = "Taxon") +
  labs(
    title = "Eukaryota",
    x = "Date",
    y = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave(
  filename = "Rel_Abundance_Euk15.png", 
  plot = last_plot(),  
  device = "png",  
  width = 10, 
  height = 7, 
  units = "in", 
  dpi = 600,  
  bg = "white" )


library(svglite)
ggsave(
  filename = "Rel_Abundance_Euk15.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)


############################
### Statistical Analyses ###
############################

## Armorphorea
Armo_data <- data_long %>%
  filter(Taxon == "Armophorea_X")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Armo_data)


## Craomonadida
Cryo_data <- data_long %>%
  filter(Taxon == "Cryomonadida")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Cryo_data)
aggregate(Abundance ~ Season, data = Cryo_data, FUN = mean)


## Dactylopodida
Dacty_data <- data_long %>%
  filter(Taxon == "Dactylopodida")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Dacty_data)
aggregate(Abundance ~ Season, data = Dacty_data, FUN = mean)

## FUngi_XX
Fung_data <- data_long %>%
  filter(Taxon == "Fungi_XX")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Fung_data)


## Himatismenida
Hima_data <- data_long %>%
  filter(Taxon == "Himatismenida")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Hima_data)


## Kinetoplastida
Kinet_data <- data_long %>%
  filter(Taxon == "Kinetoplastida")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Kinet_data)
aggregate(Abundance ~ Season, data = Kinet_data, FUN = mean)


## Microsporidiomycotina
Micro_data <- data_long %>%
  filter(Taxon == "Microsporidiomycotina")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Micro_data)


## Peritirichia
Peri_data <- data_long %>%
  filter(Taxon == "Peritrichia")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Peri_data)

## Saccharomycotina
Sac_data <- data_long %>%
  filter(Taxon == "Saccharomycotina")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Sac_data)

## Vannellida
Vann_data <- data_long %>%
  filter(Taxon == "Vannellida")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Vann_data)




#####################################################################################################
########################################## Prokaryotes ##############################################
#####################################################################################################


####################################
### Load Required Libraries
#####################################

library(ggplot2)
library(reshape2)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)


#####################################
### Step 1: Load rarefied counttables & environmental tables
#####################################

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
# select only genera, which occur in data_rarefied_RNA
prokaryote_genera_present <- unique(count_data_silva$Genus[count_data_silva$Genus %in% colnames(data_rarefied_RNA)])

# Filter count_data_silva, only for relevant orders
genus_to_order_silva_filtered <- count_data_silva %>%
  filter(Genus %in% prokaryote_genera_present) %>%
  select(Genus, Order) %>%
  distinct()

# extract relevant genera from RNA table
data_prokaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(prokaryote_genera_present))

# rename genera columns to specific order
genus_to_order_map <- setNames(genus_to_order_silva_filtered$Order, genus_to_order_silva_filtered$Genus)

# rename column
data_prokaryote_order <- data_prokaryote %>%
  pivot_longer(-SRR_ID, names_to = "Genus", values_to = "Abundance") %>%
  left_join(genus_to_order_silva_filtered, by = "Genus") %>%
  group_by(SRR_ID, Order) %>%  
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Order, values_from = Abundance, values_fill = 0)


#####################################
### Step 4: Extract Eukaryotes Genera & map at Order
#####################################
# select only genera, which occur in data_rarefied_RNA 
eukaryote_genera_present <- unique(count_data_pr2$Genus[count_data_pr2$Genus %in% colnames(data_rarefied_RNA)])

# Filter count_data_pr2, only for relevant orders
genus_to_order_pr2_filtered <- count_data_pr2 %>%
  filter(Genus %in% eukaryote_genera_present) %>%
  select(Genus, Order) %>%
  distinct()

# extrac relevant genera from RNA table
data_eukaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(eukaryote_genera_present))

# rename genera columns to specific order
genus_to_order_map <- setNames(genus_to_order_pr2_filtered$Order, genus_to_order_pr2_filtered$Genus)

# ️rename columsn
data_eukaryote_order <- data_eukaryote %>%
  pivot_longer(-SRR_ID, names_to = "Genus", values_to = "Abundance") %>%
  left_join(genus_to_order_pr2_filtered, by = "Genus") %>%
  group_by(SRR_ID, Order) %>%  
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Order, values_from = Abundance, values_fill = 0)  




#####################################
### Step 5: Remove Outlier from Metadata
#####################################

# Remove outlier samples (are already exluded from rarefied data table, so only exclude them from Env. table)
removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)

removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
data_rarefied_refseq <- data_rarefied_refseq %>% filter(!SRR_ID %in% removed_refseq_samples)
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)





#####################################
### Step 6: Set Row Names Refseq & environmental table
#####################################

data_prokaryote_order <- as.data.frame(data_prokaryote_order)
rownames(data_prokaryote_order) <- data_prokaryote_order$SRR_ID
data_prokaryote_order <- data_prokaryote_order[, -1]


#####################################
### Step 7: Relativize data
#####################################
# Compute relative abundance only for taxonomic counts
data_prokaryote_order <- data_prokaryote_order / rowSums(data_prokaryote_order)




#####################################
### Step 8: Convert to long format
#####################################
# Convert data to long format for ggplot2
data_long <- data_prokaryote_order %>%
  tibble::rownames_to_column(var = "SRR_ID") %>%
  pivot_longer(-SRR_ID, names_to = "Taxon", values_to = "Abundance")  

# Summarize top 15 most abundant taxa across all samples
top_taxa <- colSums(data_prokaryote_order) %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  head(15)

print(top_taxa)

# Group non-top taxa into "Other"
data_long <- data_long %>%
  mutate(Taxon = ifelse(Taxon %in% top_taxa, Taxon, "Other"))

# Merge with environmental metadata 
data_long <- left_join(data_long, Environmental_Data_RNA %>% select(SRR_ID, Date, Season), by = "SRR_ID")

# Konvertiere Date ins richtige Format
data_long$Date <- as.Date(data_long$Date, format = "%Y-%m-%d")




#####################################
### Step 9: Plot relative abundances over Seasons
#####################################
data_long$Taxon <- factor(data_long$Taxon, levels = c("Other", setdiff(unique(data_long$Taxon), "Other")))

# set order of seasons
data_long$Season <- factor(data_long$Season, levels = c("winter", "spring", "summer", "autumn"))


custom_colors <- c(
  "Enterobacterales" = "#606c38",  
  "Nitrososphaerales" = "#40916c",  
  "Methanosarciniales" = "#80ed99",  
  "Bacillales" = "#b7e4c7",  
  "Lactobacillales" = "#80ffdb",  
  "Burkholderiales" = "#56cfe1",  
  "Chitinophagales" = "#5390d9",  
  "Microtrichales" = "#6930c3",  
  "Bacteroidales" = "#3c096c",  
  "Pseudomonadales" = "#9d4edd", 
  "Flavobacteriales" = "#ff8500", 
  "Methanomicrobiales" = "#ff8fa3", 
  "Methanobacteriales" = "#ff4d6d", 
  "Leptospirales" = "#a4133c", 
  "Campylobacterales" = "#590d22",
  "Sphingobacteriales" = "#ff4800",
  "Woesearchaeales" = "#ff8500",
  "Lachnospirales" = "#ffaa00",
  "Absconditabacteriales SR1" = "#ffd000",
  "Pirellulales" = "#ffea00",
  "Planctomycetales" = "#e7bc91",
  "Micrococcales" = "#bc8a5f",
  "Fusobacteriales" = "#6f4518",
  "Ardenticatenales" = "#748cab",
  "GracilibacteriaX" = "#3e5c76",
  "Propionibacteriales" = "#1d2d44",
  "Candidatus Nomurabacteria" = "#a5be00",
  "Staphylococcales" = "#9d8189",
  "Saccharimonadales" = "#706677",
  "Oscillospirales" = "#a53860",
  "Other" = "#403d39"    
)



#####################################
### Step 10: Plot relative abundances over all Dates
#####################################
# Convert `Date` into correct format and sort date
data_long$Date <- as.Date(data_long$Date, format = "%Y-%m-%d")
data_long <- data_long %>% arrange(Date)  


# Find first date in month
first_dates <- data_long %>%
  group_by(format(Date, "%Y-%m-%d")) %>%
  summarise(first_date = min(Date)) %>%
  pull(first_date)

# create plot
ggplot(data_long, aes(x = as.factor(Date), y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_x_discrete(labels = format(as.Date(unique(data_long$Date)), "%Y-%m-%d")) +  
  scale_fill_manual(values = custom_colors, name = "Taxon") +
  labs(
    title = "Prokaryota",
    x = "Date",
    y = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold")
  )

ggsave(
  filename = "Rel_Abundance_Prok15.png",  
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
  filename = "Rel_Abundance_Prok15.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)


#####################################
###### Statistical testing
#####################################

####different between seasons
## Bacillales
Baci_data <- data_long %>%
  filter(Taxon == "Bacillales")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Baci_data)
aggregate(Abundance ~ Season, data = Baci_data, FUN = mean)

## Bacteroidales
Bact_data <- data_long %>%
  filter(Taxon == "Bacteroidales")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Bact_data)


## Burkholderiales
# Nur Daten für Bacillales auswählen
Burk_data <- data_long %>%
  filter(Taxon == "Burkholderiales")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Burk_data)

## Campylobacterales
Camp_data <- data_long %>%
  filter(Taxon == "Campylobacterales")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Camp_data)

## CHitinophagales
Chit_data <- data_long %>%
  filter(Taxon == "Chitinophagales")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Chit_data)

## Enterobacterales
Ent_data <- data_long %>%
  filter(Taxon == "Enterobacterales")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Ent_data)
aggregate(Abundance ~ Season, data = Ent_data, FUN = mean)

## Flavobacteriale
Flav_data <- data_long %>%
  filter(Taxon == "Flavobacteriales")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Flav_data)

## Lactobacillales
Lact_data <- data_long %>%
  filter(Taxon == "Lactobacillales")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Lact_data)
aggregate(Abundance ~ Season, data = Lact_data, FUN = mean)

## Leptospirales
Lept_data <- data_long %>%
  filter(Taxon == "Leptospirales")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Lept_data)

## Methanobacteriales
Meth_data <- data_long %>%
  filter(Taxon == "Methanobacteriales")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Meth_data)


## Methanmicrobiales
Metha_data <- data_long %>%
  filter(Taxon == "Methanomicrobiales")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Metha_data)


## Methanosarcinales
Meths_data <- data_long %>%
  filter(Taxon == "Methanosarciniales")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Meths_data)

## Microtrichales
Micro_data <- data_long %>%
  filter(Taxon == "Microtrichales")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Micro_data)

## Nitrososphaerales
Nitro_data <- data_long %>%
  filter(Taxon == "Nitrososphaerales")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Nitro_data)

## Pseudomonadales
Pseud_data <- data_long %>%
  filter(Taxon == "Pseudomonadales")

# Kruskal-Wallis-Test: Differences of abundances of Bacillales between Seasons?
kruskal.test(Abundance ~ Season, data = Pseud_data)





#####################################################################################################
########################################## ARGs #####################################################
#####################################################################################################


####################################
### Load Required Libraries
#####################################
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(tibble)

#####################################
### Step 1: Load rarefied counttables & environmental tables
#####################################

# ARG from 16S normalized table
ARG_data <- read_csv("Master/Metatranscriptomics/Excel_lists/ARG/ARG_data_combined.csv")

Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")


#####################################
### Step 2: Prepare ARG table
#####################################

# Convert into Long-Format
ARG_long <- ARG_data %>%
  pivot_longer(cols = -Family, names_to = "SRR_ID", values_to = "Abundance")

# Convert back to Wide-Format, but with SRR_ID as rowname
ARG_transformed <- ARG_long %>%
  pivot_wider(names_from = Family, values_from = Abundance, values_fill = 0)

# set SRR_ID as rowname
ARG_matrix <- ARG_transformed %>%
  column_to_rownames("SRR_ID")


#####################################
### Step 3: Remove outliers
#####################################

removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)
ARG_matrix <- ARG_matrix[!rownames(ARG_matrix) %in% removed_refseq_samples, ]


#####################################
### Step 4: Relativize data (per Sample)
#####################################
ARG_matrix <- ARG_matrix / rowSums(ARG_matrix)


#####################################
### Step 5: Convert to long format
#####################################
# convert into Long-Format for ggplot2
data_long <- ARG_matrix %>%
  tibble::rownames_to_column(var = "SRR_ID") %>%
  pivot_longer(-SRR_ID, names_to = "Taxon", values_to = "Abundance")

# extract 10 most abundant ARG families
top_taxa <- colSums(ARG_matrix) %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  head(10)

print(top_taxa)  

# group not abundant into "Other"
data_long <- data_long %>%
  mutate(Taxon = ifelse(Taxon %in% top_taxa, Taxon, "Other"))

# show taxa which are in Other
not_top10_taxa <- setdiff(colnames(ARG_matrix), top_taxa)
print(not_top10_taxa)


# Merge with environmental data
data_long <- left_join(data_long, Environmental_Data_DNA %>% select(SRR_ID, Date, Season), by = "SRR_ID")

# date as correct date
data_long <- data_long %>%
  mutate(Date = as.Date(Date, format = "%Y-%m-%d")) %>%
  arrange(Date, SRR_ID)  


#####################################
### Step 6: Plot relative abundances over Seasons
#####################################

# set order of seasons
data_long$Season <- factor(data_long$Season, levels = c("winter", "spring", "summer", "autumn"))


custom_colors <- c(
  "Aminoglycoside" = "#80ffdb",  
  #"antibacterial_fatty_acid" = "",  
  "Bacitracin" = "#56cfe1",  
  "Beta_lactam" = "#5390d9",  
  #"bicyclomycin" = "",  
  #"bleomycin" = "",  
  #"chloramphenicol" = "",  
  #"defensin" = "",  
  #"florfenicol" = "", 
  # "fosfomycin" = "", 
  "Macrolide-lincosamide-streptogramin" = "#6930c3", 
  "Multidrug" = "#9d4edd", 
  "Mupirocin" = "#e7c6ff", 
  #"novobiocin" = "",
  "Polymyxin" = "#ff4800",
  #"quinolone" = "",
  "Rifamycin" = "#ff8500",
  "Sulfonamide" = "#e7bc91",
  "Tetracycline" = "#a53860",
  #"trimethoprim" = "#e7bc91",
  "Other" = "#403d39"
)

# Set other at last position
data_long$Taxon <- factor(data_long$Taxon, levels = c("Other", setdiff(unique(data_long$Taxon), "Other")))


#####################################
### Step 7: Plot relative abundances over all Dates
#####################################

ggplot(data_long, aes(x = as.factor(Date), y = Abundance, fill = Taxon)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) + 
  scale_x_discrete(labels = format(as.Date(unique(data_long$Date)), "%Y-%m-%d")) +  
  scale_fill_manual(values = custom_colors, name = "Taxon") +
  labs(
    title = "ARG",
    x = "Date",
    y = "Relative Abundance"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
    legend.position = "bottom",
    plot.title = element_text(size = 14, face = "bold")
  )


ggsave(
  filename = "Rel_Abundance_ARG10.png",
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
  filename = "Rel_Abundance_ARG10.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)


#####################################
###### Statistical testing
#####################################
# Normal distribution?
shapiro.test(data_long$Abundance) 

#### not different between seasons
## Aminoglycoside
Amino_data <- data_long %>%
  filter(Taxon == "Aminoglycoside")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Amino_data)

## Bacitracin
Baci_data <- data_long %>%
  filter(Taxon == "Bacitracin")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Baci_data)
aggregate(Abundance ~ Season, data = Baci_data, FUN = mean)

## Beta-lactam
Beta_data <- data_long %>%
  filter(Taxon == "Beta_lactam")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Beta_data)
aggregate(Abundance ~ Season, data = Beta_data, FUN = mean)

##Macrolide
Macro_data <- data_long %>%
  filter(Taxon == "Macrolide-lincosamide-streptogramin")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Macro_data)

##Multidrug
Multi_data <- data_long %>%
  filter(Taxon == "Multidrug")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Multi_data)
aggregate(Abundance ~ Season, data = Multi_data, FUN = mean)

##Mupirocin
Mui_data <- data_long %>%
  filter(Taxon == "Mupirocin")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Mui_data)


##Polymyxin
Poly_data <- data_long %>%
  filter(Taxon == "Polymyxin")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Poly_data)


## Rifamycin
Rifa_data <- data_long %>%
  filter(Taxon == "Rifamycin")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Rifa_data)



## Sulfonamide
Sulfo_data <- data_long %>%
  filter(Taxon == "Sulfonamide")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Sulfo_data)


##Tetracycline
Tetra_data <- data_long %>%
  filter(Taxon == "Tetracycline")

# Kruskal-Wallis-Test: Differences of abundances of Kyanoviridae between Seasons?
kruskal.test(Abundance ~ Season, data = Tetra_data)
