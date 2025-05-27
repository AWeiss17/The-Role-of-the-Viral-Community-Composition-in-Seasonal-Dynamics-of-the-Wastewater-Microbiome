###########################################################################################################
############################## Spearman Correlation Viruses - ARGs ########################################
###########################################################################################################


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

# Table ARG & DNA rarefied
ARG_data <- read_csv("Master/Metatranscriptomics/Excel_lists/ARG/ARG_data_combined.csv")
data_rarefied_refseq <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Rarefied_count_data_refseq_486.csv")

# Environmental RNA table
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")



#####################################
### Step 2: Prepare tables
#####################################
data_rarefied_refseq <- data_rarefied_refseq %>% rename(SRR_ID = ...1)

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
removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)
ARG_transformed <- ARG_transformed %>% filter(!SRR_ID %in% removed_refseq_samples)
data_rarefied_refseq <- data_rarefied_refseq %>% filter(!SRR_ID %in% removed_refseq_samples)

# select only neccessary Metadata columns
Environmental_Data_DNA <- Environmental_Data_DNA %>%
  select(SRR_ID, Date, Season, Temp_manual, pH_online, Oxygen_manual_mg_L, Conductivity_manual_mS_cm, Temp_Air, Oxygen_sat_manual, NO3_online_mg.L, PO4_online_mg.L)  # Nur relevante Umweltvariablen behalten


#####################################
### Step 3: Format Data
#####################################

# Set Data as rowname
data_rarefied_refseq <- column_to_rownames(data_rarefied_refseq, var = "SRR_ID")
ARG_matrix <- column_to_rownames(ARG_transformed, var = "SRR_ID")

# relativize ARG data
ARG_matrix <- ARG_matrix / rowSums(ARG_matrix)
#data_rarefied_refseq <- data_rarefied_refseq / rowSums(data_rarefied_refseq)

# transform back to data frame
ARG_matrix <- as.data.frame(ARG_matrix)
ARG_matrix <- ARG_matrix[match(rownames(data_rarefied_refseq), rownames(ARG_matrix)), ]
identical(rownames(data_rarefied_refseq), rownames(ARG_matrix)) 



#####################################
### Step 4: select top 10 taxa from Virus and Prokaryotes
#####################################

# Claculate sum of every viral family across samples
top_taxa_Vir <- colSums(data_rarefied_refseq) %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  head(10)  


top_taxa_arg <- colSums(ARG_matrix) %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  head(10)  


# Filtere die Data
data_rarefied_top10 <- data_rarefied_refseq[, colnames(data_rarefied_refseq) %in% top_taxa_Vir]
data_arg_top10 <- ARG_matrix[, colnames(ARG_matrix) %in% top_taxa_arg]



#####################################
### Step 5: Calculate Spearman Correlation 
#####################################
# Initialisize Matrix for correlations
correlation_matrix <- matrix(NA, 
                             nrow = ncol(data_rarefied_top10), 
                             ncol = ncol(data_arg_top10),
                             dimnames = list(colnames(data_rarefied_top10), colnames(data_arg_top10)))


p_value_matrix <- correlation_matrix  

# Calculate Spearman-Correlations + p-values
for (virus in colnames(data_rarefied_top10)) {
  for (arg in colnames(data_arg_top10)) {
    test <- cor.test(data_rarefied_top10[[virus]], data_arg_top10[[arg]], method = "spearman")
    correlation_matrix[virus, arg] <- test$estimate 
    p_value_matrix[virus, arg] <- test$p.value  
  }
}



#####################################
### Step 6: Plot Heatmap 
#####################################
# Correlations und p-values
melted_correlation <- melt(correlation_matrix, na.rm = TRUE)
melted_p_values <- melt(p_value_matrix, na.rm = TRUE)

melted_data <- left_join(melted_correlation, melted_p_values, by = c("Var1", "Var2"), suffix = c("_cor", "_pval"))

# sign significant values with asterrisk
melted_data <- melted_data %>%
  mutate(significance = case_when(
    value_pval < 0.001 ~ "***",  
    value_pval < 0.01  ~ "**",   
    value_pval < 0.05  ~ "*",    
    TRUE ~ ""                    
  ))



# calculate min/max
fill_limits <- range(melted_data$value_cor, na.rm = TRUE)

# Heatmap 
heatmap_plot <- ggplot(melted_data, aes(x = Var2, y = Var1, fill = value_cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#d90429", mid = "white", high = "#0096c7", 
                       midpoint = 0, limits = fill_limits, name = "Correlation", na.value = "white") +
  geom_text(aes(label = significance), color = "black", size = 5) +  # Signifikante Werte mit Stern markieren
  labs(title = "Spearman Correlation",
       x = "ARG Families",
       y = "Virus Families") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(heatmap_plot)

ggsave("Virus_ARG_Top10_Heatmap.pdf", plot = heatmap_plot, width = 8, height = 10)

library(writexl)
# only significant Correlations (p < 0.05)
sig_only <- melted_data %>%
  filter(value_pval < 0.05)


write_xlsx(sig_only, "Virus_ARG_Significant_Correlations.xlsx")


##############################
##### Shapiro-Wilk
##############################
# Test for normality

normality_summary <- function(df) {
  apply(df, 2, function(x) shapiro.test(x)$p.value)
}

normality_summary(data_rarefied_top10)
normality_summary(data_arg_top10)


ggsave(
  filename = "Spearman_ARG.png",  
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
  filename = "Spearman_ARG.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)





###############################################################################################
################### Spearman Correlation Viruses - Prokaryotes ################################
###############################################################################################


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
# select only genera auswählen in data_rarefied_RNA
prokaryote_genera_present <- unique(count_data_silva$Genus[count_data_silva$Genus %in% colnames(data_rarefied_RNA)])

genus_to_order_silva_filtered <- count_data_silva %>%
  filter(Genus %in% prokaryote_genera_present) %>%
  select(Genus, Order) %>%
  distinct()

data_prokaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(prokaryote_genera_present))

# rename genera columns
genus_to_order_map <- setNames(genus_to_order_silva_filtered$Order, genus_to_order_silva_filtered$Genus)

# rename columns
data_prokaryote_order <- data_prokaryote %>%
  pivot_longer(-SRR_ID, names_to = "Genus", values_to = "Abundance") %>%
  left_join(genus_to_order_silva_filtered, by = "Genus") %>%
  group_by(SRR_ID, Order) %>% 
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Order, values_from = Abundance, values_fill = 0)  

# Set SRR_ID as rowname and convert to Matrix
Prok_data_matrix <- data_prokaryote_order %>%
  column_to_rownames(var = "SRR_ID") %>%
  as.matrix()

# Normalization of the data is neccessary!
# Although we use the rarefied table, which is already "normalized" 
# BUT! We divide pr2 and silva dataset! So new relativization is neccessary!
# Relativize per sample
Prok_data_rel <- Prok_data_matrix / rowSums(Prok_data_matrix)

# Matrix back in DataFrame
Prok_data_table <- as.data.frame(Prok_data_rel)

# SRR_ID rownames was column
Prok_data_table <- Prok_data_table %>%
  rownames_to_column(var = "SRR_ID")



#####################################
### Step 4: Extract Eukaryotes Genera & map at Order
#####################################
# only genera in data_rarefied_RNA 
eukaryote_genera_present <- unique(count_data_pr2$Genus[count_data_pr2$Genus %in% colnames(data_rarefied_RNA)])

# Filtere count_data_pr2
genus_to_order_pr2_filtered <- count_data_pr2 %>%
  filter(Genus %in% eukaryote_genera_present) %>%
  select(Genus, Order) %>%
  distinct()

data_eukaryote <- data_rarefied_RNA %>%
  select(SRR_ID, any_of(eukaryote_genera_present))

genus_to_order_map <- setNames(genus_to_order_pr2_filtered$Order, genus_to_order_pr2_filtered$Genus)

# ️rename column
data_eukaryote_order <- data_eukaryote %>%
  pivot_longer(-SRR_ID, names_to = "Genus", values_to = "Abundance") %>%
  left_join(genus_to_order_pr2_filtered, by = "Genus") %>%
  group_by(SRR_ID, Order) %>%  
  summarise(Abundance = sum(Abundance, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Order, values_from = Abundance, values_fill = 0)  

# Set SRR_ID as rowname and convert to Matrix
Euk_data_matrix <- data_eukaryote_order %>%
  column_to_rownames(var = "SRR_ID") %>%
  as.matrix()

# Normalization of the data is neccessary!
# Although we use the rarefied table, which is already "normalized" 
# BUT! We divide pr2 and silva dataset! So new relativization is neccessary!
# Relativize per sample (row)
Euk_data_rel <- Euk_data_matrix / rowSums(Euk_data_matrix)

# Matrix back in DataFrame
Euk_data_table <- as.data.frame(Euk_data_rel)

# SRR_ID from rownames as column
Euk_data_table <- Euk_data_table %>%
  rownames_to_column(var = "SRR_ID")


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
### Step 6: Set Date as rowname (for later combination of tables)
#####################################

data_rarefied_refseq <- data_rarefied_refseq %>%
  left_join(Environmental_Data_DNA %>% select(SRR_ID, Date), by = "SRR_ID") %>%
  select(Date, everything(), -SRR_ID)  


Prok_data_table <- Prok_data_table %>%
  left_join(Environmental_Data_RNA %>% select(SRR_ID, Date), by = "SRR_ID") %>%
  select(Date, everything(), -SRR_ID)  

Euk_data_table <- Euk_data_table %>%
  left_join(Environmental_Data_RNA %>% select(SRR_ID, Date), by = "SRR_ID") %>%
  select(Date, everything(), -SRR_ID)  




#####################################
### Step 7: Format Data
#####################################

# Set Date rowname
data_rarefied_refseq <- column_to_rownames(data_rarefied_refseq, var = "Date")
Prok_data_table <- column_to_rownames(Prok_data_table, var = "Date")
Euk_data_table <- column_to_rownames(Euk_data_table, var = "Date")


# check similar order
Prok_data_table <- Prok_data_table[match(rownames(data_rarefied_refseq), rownames(Prok_data_table)), ]
Euk_data_table <- Euk_data_table[match(rownames(data_rarefied_refseq), rownames(Euk_data_table)), ]

identical(rownames(data_rarefied_refseq), rownames(Prok_data_table)) 
identical(rownames(data_rarefied_refseq), rownames(Euk_data_table)) 



####################################
### Step 8: select top 10 taxa from Virus and Prokaryotes
#####################################

# Calculate sum of varus families across all samples
top_taxa_Vir <- colSums(data_rarefied_refseq) %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  head(10)  


top_taxa_prok <- colSums(Prok_data_table) %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  head(15)  


top_taxa_euk <- colSums(Euk_data_table) %>%
  sort(decreasing = TRUE) %>%
  names() %>%
  head(15) 



# Filter data top Top 10 Taxa
data_rarefied_top10 <- data_rarefied_refseq[, colnames(data_rarefied_refseq) %in% top_taxa_Vir]

data_prokaryote_top15 <- Prok_data_table[, colnames(Prok_data_table) %in% top_taxa_prok]




#####################################
### Step 9: Calculate Spearman Correlation 
#####################################
# Initialize Matrix for Correlationen
correlation_matrix <- matrix(NA, 
                             nrow = ncol(data_rarefied_top10), 
                             ncol = ncol(data_prokaryote_top15),
                             dimnames = list(colnames(data_rarefied_top10), colnames(data_prokaryote_top15)))

p_value_matrix <- correlation_matrix  

# Calculate Spearman-Correlationen + p-values
for (virus in colnames(data_rarefied_top10)) {
  for (pro in colnames(data_prokaryote_top15)) {
    test <- cor.test(data_rarefied_top10[[virus]], data_prokaryote_top15[[pro]], method = "spearman")
    correlation_matrix[virus, pro] <- test$estimate  
    p_value_matrix[virus, pro] <- test$p.value  
  }
}



#####################################
### Step 6: Plot Heatmap 
#####################################
# Melt correlations and p-values
melted_correlation <- melt(correlation_matrix, na.rm = TRUE)
melted_p_values <- melt(p_value_matrix, na.rm = TRUE)

# add p-values
melted_data <- left_join(melted_correlation, melted_p_values, by = c("Var1", "Var2"), suffix = c("_cor", "_pval"))

# sign significant correlations with asterrisk
melted_data <- melted_data %>%
  mutate(significance = case_when(
    value_pval < 0.001 ~ "***",  
    value_pval < 0.01  ~ "**",   
    value_pval < 0.05  ~ "*",    
    TRUE ~ ""                    
  ))



# Calculate Min-/Max-values
fill_limits <- range(melted_data$value_cor, na.rm = TRUE)

# Heatmap
heatmap_plot <- ggplot(melted_data, aes(x = Var2, y = Var1, fill = value_cor)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "#d90429", mid = "white", high = "#0096c7", 
                       midpoint = 0, limits = fill_limits, name = "Correlation", na.value = "white") +
  geom_text(aes(label = significance), color = "black", size = 5) + 
  labs(title = "Spearman Correlation: Virus & Prok.",
       x = "Prokaryote Orders",
       y = "Virus Families") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(heatmap_plot)


# save plot
ggsave("Virus_Prok_order_Heatmap.pdf", plot = heatmap_plot, width = 8, height = 10)


library(writexl)
# only significant values(p < 0.05)
sig_only <- melted_data %>%
  filter(value_pval < 0.05)

# Save in Excel
write_xlsx(sig_only, "Virus_Prok_Significant_Correlations.xlsx")

##############################
##### Shapiro-Wilk
##############################
# Test for normality

normality_summary <- function(df) {
  apply(df, 2, function(x) shapiro.test(x)$p.value)
}

normality_summary(data_rarefied_top10)
normality_summary(data_prokaryote_top15)


ggsave(
  filename = "Spearman_Prok.png",  
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
  filename = "Spearman_Prok.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)

