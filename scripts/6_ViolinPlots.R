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
library(lubridate)
library(car)
library(svglite)

# MG and MT complete reads
data <- read_delim("Master/Metatranscriptomics/Excel_lists/accession_with_some_info.csv", 
                   delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Load Environmental Data (DNA)
Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")


# Load Viral count data (not rarefied)
data_refseq_wide_SRR <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/count_data_refseq_wide_SRR.csv")

# Remove unnecessary columns
data_refseq_wide_SRR <- data_refseq_wide_SRR %>%
  rename(SRR_ID = ...1)  


#########################
##### Prepare data table
#########################
# remove samples
removed_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506", "SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
data <- data %>% filter(!SRR_ID %in% removed_samples)

# remove samples from Env data DNA
removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)
data_refseq_wide_SRR <- data_refseq_wide_SRR %>% filter(!SRR_ID %in% removed_refseq_samples)

Environmental_Data_DNA <- Environmental_Data_DNA %>%
  select(SRR_ID, Date, Season)  # Keep only necessary environmental variables

# Ensure Date is in Date format
data <- data %>%
  mutate(Date = dmy(Date))


###################
#### Metagenomic
###################
# Filter only Metagenomic data
data_MG <- data %>% filter(LibrarySource == "METAGENOMIC")

# to see ho much Metagenomic reads in total were obtained across the selected samples
total_metagenomic_reads <- sum(data_MG$Reads, na.rm = TRUE)

total_vir_reads <- sum(data_refseq_wide_SRR[,-1], na.rm = TRUE)

vir_percentage <- (total_vir_reads / total_metagenomic_reads) * 100
# combine tables MG with Env
data_MG <- merge(
  data_MG,
  Environmental_Data_DNA,
  by = "SRR_ID",
  all.x = TRUE
)

# delet one Date column
data_MG <- data_MG %>%
  select(-Date.y) %>%
  rename(Date = Date.x)



########################
#### Prepare Virus table
########################
# Sum sequences per sample
virus_counts <- data_refseq_wide_SRR %>%
  mutate(Sequences = rowSums(select(., -SRR_ID))) %>%  # Sum all viral sequences per sample
  select(SRR_ID, Sequences)  # Keep only SRR_ID and summed sequences

# combine virus table wit MG
virus_data_rel<- left_join(
  data_MG,
  virus_counts,
  by = "SRR_ID"
)

# calculae relative amount in MG dataset
virus_data_rel <- virus_data_rel %>%
  mutate(rel_virus_reads = Sequences / Reads)

# in percent
virus_data_rel <- virus_data_rel %>%
  mutate(rel_virus_reads_percent = rel_virus_reads * 100)


# Convert Season to an ordered factor
virus_data_rel$Season <- factor(virus_data_rel$Season, levels = c("spring", "summer", "autumn", "winter"))

###Violin-Plot per Season
ggplot(virus_data_rel, aes(x = Season, y = rel_virus_reads_percent, fill = Season)) +
  geom_violin(alpha = 0.5) +
  geom_jitter(shape = 16, position = position_jitter(0.2), aes(color = Season), size = 1.5, alpha = 0.6) +
  labs(
    title = "Viruses",
    x = "",
    y = "Viral content [%]",
    fill = "Season",
    color = "Season"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  
  scale_fill_manual(values = c("spring" = "#62b6cb", "summer" = "#ff7d00", "autumn" = "#8c1c13", "winter" = "#134074")) +
  scale_color_manual(values = c("spring" = "#62b6cb", "summer" = "#ff7d00", "autumn" = "#8c1c13", "winter" = "#134074")) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12, face = "bold")
  )

ggsave(
  filename = "Violine_Vir.png",  
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
  filename = "Violine_Vir.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)

##################################
####### Step 2: statistic Tests
##################################

shapiro.test(virus_data_rel$rel_virus_reads_percent)  # Test auf Normalverteilung
kruskal.test(rel_virus_reads_percent ~ Season, data = virus_data_rel)





#####################################################################################################
########################################## Prokaryotes ##############################################
#####################################################################################################


####################################
### Load Required Libraries
#####################################

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(car)
library(svglite)

# MG and MT complete reads
data <- read_delim("Master/Metatranscriptomics/Excel_lists/accession_with_some_info.csv", 
                   delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Load Prokaryotic count data (RNA) (not rarefied)
data_silva_wide_SRR <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/count_data_silva_wide_SRR.csv")

# Load Environmental Data (RNA)
Environmental_Data_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metatranscriptomics.csv")


# Remove unnecessary columns
data_silva_wide_SRR <- data_silva_wide_SRR %>%
  rename(SRR_ID = ...1)  

#########################
##### Prepare data table
#########################
# remove samples
removed_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506", "SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
data <- data %>% filter(!SRR_ID %in% removed_samples)

Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season)  # Keep only necessary environmental variables


# Remove outlier samples
removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
data_silva_wide_SRR <- data_silva_wide_SRR %>% filter(!SRR_ID %in% removed_rna_samples)
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)

# Ensure Date is in Date format
data <- data %>%
  mutate(Date = dmy(Date))



###################
#### Metatranscriptomic
###################
data_MT <- data %>% filter(LibrarySource == "METATRANSCRIPTOMIC")

# Calculate proportion of Prokaryotic reads in total
# to see ho much Metatranscriptomic reads in total were obtained across the selected samples
total_metatrans_reads <- sum(data_MT$Reads, na.rm = TRUE)

total_prok_reads <- sum(data_silva_wide_SRR[,-1], na.rm = TRUE)

prok_percentage <- (total_prok_reads / total_metatrans_reads) * 100


# combine tables MG with Env
data_MT <- merge(
  data_MT,
  Environmental_Data_RNA,
  by = "SRR_ID",
  all.x = TRUE
)

# delet one Date column
data_MT <- data_MT %>%
  select(-Date.y) %>%
  rename(Date = Date.x)



########################
#### Prepare Prokaryote table
########################
# Sum sequences per sample
prokaryote_counts <- data_silva_wide_SRR %>%
  mutate(Sequences = rowSums(select(., -SRR_ID))) %>%  
  select(SRR_ID, Sequences) 


# combine virus table wit MG
prokaryote_counts_rel<- left_join(
  data_MT,
  prokaryote_counts,
  by = "SRR_ID"
)

# calculae relative amount in MG dataset
prokaryote_counts_rel <- prokaryote_counts_rel %>%
  mutate(rel_prok_reads = Sequences / Reads)

# in percent
prokaryote_counts_rel <- prokaryote_counts_rel %>%
  mutate(rel_prok_reads_percent = rel_prok_reads * 100)


# Convert Season to an ordered factor
prokaryote_counts_rel$Season <- factor(prokaryote_counts_rel$Season, levels = c("spring", "summer", "autumn", "winter"))


###Violin-Plot per Season
ggplot(prokaryote_counts_rel, aes(x = Season, y = rel_prok_reads_percent, fill = Season)) +
  geom_violin(alpha = 0.5) +
  geom_jitter(shape = 16, position = position_jitter(0.2), aes(color = Season), size = 1.5, alpha = 0.6) +
  labs(
    title = "Virus",
    x = "",
    y = "Viral content [%]",
    fill = "Season",
    color = "Season"
  ) +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +  # Optional für %-Format
  scale_fill_manual(values = c("spring" = "#62b6cb", "summer" = "#ff7d00", "autumn" = "#8c1c13", "winter" = "#134074")) +
  scale_color_manual(values = c("spring" = "#62b6cb", "summer" = "#ff7d00", "autumn" = "#8c1c13", "winter" = "#134074")) +
  theme_minimal() +
  theme(
    legend.title = element_text(size = 12, face = "bold")
  )

ggsave(
  filename = "Violine_prok.png",  
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
  filename = "Violine_prok.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)


##################################
####### Step 2: statistic Tests
##################################

shapiro.test(prokaryote_counts_rel$rel_prok_reads_percent)  # Test auf Normalverteilung
kruskal.test(rel_prok_reads_percent~ Season, data = prokaryote_counts_rel)







#####################################################################################################
########################################## Eukaryotes ###############################################
#####################################################################################################


####################################
### Load Required Libraries
#####################################
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(lubridate)
library(car)
library(svglite)

# MG and MT complete reads
data <- read_delim("Master/Metatranscriptomics/Excel_lists/accession_with_some_info.csv", 
                   delim = ";", escape_double = FALSE, trim_ws = TRUE)

# Load Eukaryotic count data (RNA) (not rarefied)
data_pr2_wide_SRR <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/count_data_pr2_wide_SRR.csv")

# Load Environmental Data (RNA)
Environmental_Data_RNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metatranscriptomics.csv")


# Remove unnecessary columns
data_pr2_wide_SRR <- data_pr2_wide_SRR %>%
  rename(SRR_ID = ...1)  # Set the first column as SRR_ID

#########################
##### Prepare data table
#########################
# remove samples
removed_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506", "SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
data <- data %>% filter(!SRR_ID %in% removed_samples)

Environmental_Data_RNA <- Environmental_Data_RNA %>%
  select(SRR_ID, Date, Season)  # Keep only necessary environmental variables


# Remove outlier samples
removed_rna_samples <- c("SRR9006565","SRR9006522","SRR9006571", "SRR9006493", "SRR9006496", "SRR9006534")
data_pr2_wide_SRR <- data_pr2_wide_SRR %>% filter(!SRR_ID %in% removed_rna_samples)
Environmental_Data_RNA <- Environmental_Data_RNA %>% filter(!SRR_ID %in% removed_rna_samples)

# Ensure Date is in Date format
data <- data %>%
  mutate(Date = dmy(Date))



###################
#### Metatranscriptomic
###################
data_MT <- data %>% filter(LibrarySource == "METATRANSCRIPTOMIC")

# Calculate proportion of Prokaryotic reads in total
# to see ho much Metatranscriptomic reads in total were obtained across the selected samples
total_metatrans_reads <- sum(data_MT$Reads, na.rm = TRUE)

total_euk_reads <- sum(data_pr2_wide_SRR[,-1], na.rm = TRUE)

euk_percentage <- (total_euk_reads / total_metatrans_reads) * 100


# combine tables MG with Env
data_MT <- merge(
  data_MT,
  Environmental_Data_RNA,
  by = "SRR_ID",
  all.x = TRUE
)

# delet one Date column
data_MT <- data_MT %>%
  select(-Date.y) %>%
  rename(Date = Date.x)



########################
#### Prepare Eukaryote table
########################
# Sum sequences per sample
eukaryote_counts <- data_pr2_wide_SRR %>%
  mutate(Sequences = rowSums(select(., -SRR_ID))) %>%  
  select(SRR_ID, Sequences) 


# combine virus table wit MG
eukaryote_counts_rel<- left_join(
  data_MT,
  eukaryote_counts,
  by = "SRR_ID"
)

# calculae relative amount in MG dataset
eukaryote_counts_rel <- eukaryote_counts_rel %>%
  mutate(rel_euk_reads = Sequences / Reads)

# in percent
eukaryote_counts_rel <- eukaryote_counts_rel %>%
  mutate(rel_euk_reads_percent = rel_euk_reads * 100)


# Convert Season to an ordered factor
eukaryote_counts_rel$Season <- factor(eukaryote_counts_rel$Season, levels = c("spring", "summer", "autumn", "winter"))


###Violin-Plot per Season
ggplot(eukaryote_counts_rel, aes(x = Season, y = rel_euk_reads_percent , fill = Season)) +
  geom_violin(alpha = 0.5) +  
  geom_jitter(shape = 16, position = position_jitter(0.2), aes(color = Season), size = 1.5, alpha = 0.6) +  
  labs(title = "Eukaryotes", 
       x = "", 
       y = "Relative abundance [%]",
       fill = "Season",  
       color = "Season") +  
  scale_fill_manual(values = c("spring" = "#62b6cb", "summer" = "#ff7d00", "autumn" = "#8c1c13", "winter" = "#134074")) +
  scale_color_manual(values = c("spring" = "#62b6cb", "summer" = "#ff7d00", "autumn" = "#8c1c13", "winter" = "#134074")) +
  theme_minimal() +
  theme(
    legend.title = element_text(size =12, face = "bold")  
  )


ggsave(
  filename = "Violine_euk.png",  
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
  filename = "Violine_euk.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)


##################################
####### Step 2: statistic Tests
##################################

shapiro.test(eukaryote_counts_rel$rel_euk_reads_percent)  # Test auf Normalverteilung
kruskal.test(rel_euk_reads_percent~ Season, data = eukaryote_counts_rel)




#####################################################################################################
########################################## ARGs #####################################################
#####################################################################################################


####################################
### Load Required Libraries
#####################################
#### without Normalization per sample!

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)


#####################################
### Step 1: Load rarefied counttables & environmental tables
#####################################

# ARG from 16S normalized table
ARG_data <- read_csv("Master/Metatranscriptomics/Excel_lists/ARG/ARG_data_combined.csv")

Environmental_Data_DNA <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/Environmental_Data_Metagenomics.csv")


#####################################
### Step 2: Prepare ARG table
#####################################
# 1. Long Format
ARG_long <- ARG_data %>%
  pivot_longer(cols = -Family, names_to = "SRR_ID", values_to = "Abundance")

# 2. Wide Format (SRR_ID bleibt als Spalte!)
ARG_table <- ARG_long %>%
  pivot_wider(names_from = Family, values_from = Abundance, values_fill = 0)

# 3. Remove outliers directly
removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558", "SRR9006580", "SRR9006587", "SRR9006506")
ARG_table <- ARG_table %>% filter(!SRR_ID %in% removed_refseq_samples)
Environmental_Data_DNA <- Environmental_Data_DNA %>% filter(!SRR_ID %in% removed_refseq_samples)

# Sum sequences per sample
ARG_counts <- ARG_table %>%
  mutate(ARG_content = rowSums(select(., -SRR_ID))) %>% 
  select(SRR_ID, ARG_content) 

# Merge with environmental data
ARG_counts <- left_join(ARG_counts, Environmental_Data_DNA, by = "SRR_ID")

# Ensure Date is in Date format
ARG_counts$Date <- as.Date(ARG_counts$Date)

# Convert Season to an ordered factor
ARG_counts$Season <- factor(ARG_counts$Season, levels = c("spring", "summer", "autumn", "winter"))

# Sort data by Date
ARG_counts <- ARG_counts %>%
  arrange(Date)

# in percent
#ARG_counts <- ARG_counts %>%
# mutate(ARG_content_percent = ARG_content * 100)


###Violin-Plot per Season
ggplot(ARG_counts, aes(x = Season, y = ARG_content, fill = Season)) +
  geom_violin(alpha = 0.5) +  
  geom_jitter(shape = 16, position = position_jitter(0.2), aes(color = Season), size = 1.5, alpha = 0.6) +  
  labs(title = "Seasonal variation in ARG content", 
       x = "", 
       y = "ARGs per 16S rRNA gene",
       fill = "Season",  
       color = "Season") +  
  scale_fill_manual(values = c("spring" = "#62b6cb", "summer" = "#ff7d00", "autumn" = "#8c1c13", "winter" = "#134074")) +
  scale_color_manual(values = c("spring" = "#62b6cb", "summer" = "#ff7d00", "autumn" = "#8c1c13", "winter" = "#134074")) +
  theme_minimal() +
  theme(
    legend.title = element_text(size =12, face = "bold")  
  )

ggsave(
  filename = "Violine_ARG.png",  
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
  filename = "Violin_ARG.svg",  
  plot = last_plot(),  
  device = "svg",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600,  
  bg = "white" 
)

##################################
####### Step 2: statistic Tests
##################################
library(car)
shapiro.test(ARG_counts$ARG_content)  
kruskal.test(ARG_content ~ Season, data = ARG_counts)
