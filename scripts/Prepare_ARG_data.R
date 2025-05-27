####################################
### Load Required Libraries
#####################################
library(dplyr)
library(readr)



#####################################
### Step 1: Load rarefied counttables & environmental tables
#####################################
# Here I combine 2 tables of ARG Data.
# I will use the Table with the genes and sum all genes of the same family.
# Then do the same with second table
# do some renamings of Family Names
# combine the tables -> call it ARG_data to easily match with all other scripts


# ARG data only 3 samples
#ARG_data_Add <- read.table("Master/Metatranscriptomics/Excel_lists/ARG/3_samples_Normalised_16_Results_dot.txt" , header = TRUE, sep = "\t", stringsAsFactors = FALSE)
#write.csv(ARG_data_Add, "Master/Metatranscriptomics/Excel_lists/ARG/3_samples_Normalised_16_Results_dot.csv", row.names = FALSE)

# ARG data only 3 samples csv table
ARG_data_Add <- read_csv("Master/Metatranscriptomics/Excel_lists/ARG/3_samples_Normalised_16_Results_dot.csv")
unique_families <- unique(ARG_data_Add$Familiy)
print(unique_families)

# ARG from 16S normalized table
Genes_normalized_16S <- read_csv("Master/Metatranscriptomics/Excel_lists/ARG/Genes_normalized_16S.csv")
unique_families <- unique(Genes_normalized_16S$Family)
print(unique_families)



#####################################
### Step 2: Sum all genes from same Family (ARG_data_Add)
#####################################
# Sum values for each family
ARG_data_Add <- ARG_data_Add %>%
  rename(Family = Familiy)

ARG_data_Add_summarized <- ARG_data_Add %>%
  group_by(Family) %>%  # Group by family (Caution mistiake! "Familiy" vs. "Family")
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>%  
  ungroup()  


# Sum values for each family
Gene_data_summarized <- Genes_normalized_16S %>%
  group_by(Family) %>%  # Group by family (Caution mistiake! "Familiy" vs. "Family")
  summarise(across(where(is.numeric), sum, na.rm = TRUE)) %>% 

# Show the different Family names
unique_families <- unique(ARG_data_Add_summarized$Family)
print(unique_families)

unique_families <- unique(Gene_data_summarized$Family)
print(unique_families)



#####################################
### Step 3: Rename Family names
#####################################
Gene_data_summarized <- Gene_data_summarized %>%
  mutate(Family = recode(Family,
                         "aminoglycoside" = "Aminoglycoside", 
                         "antibacterial_fatty_acid" = "Antibacterial_fatty_acid", 
                         "bacitracin" = "Bacitracin",
                         "beta_lactam" = "Beta_lactam",
                         "bicyclomycin" = "Bicyclomycin",
                         "bleomycin" = "Bleomycin",
                         "chloramphenicol" = "Chloramphenicol",
                         "defensin" = "Defensin",
                         "factumycin" = "Factumycin",
                         "florfenicol" = "Florfenicol",
                         "fosfomycin" = "Fosfomycin",
                         "macrolide-lincosamide-streptogramin" = "Macrolide-lincosamide-streptogramin", 
                         "multidrug" = "Multidrug", 
                         "mupirocin" = "Mupirocin",
                         "novobiocin" = "Novobiocin", 
                         "other peptidic antibiotic" = "other peptidic antibiotic",
                         "pleuromutilin_tiamulin" = "Pleuromutilin_tiamulin",
                         "polymyxin" = "Polymyxin",
                         "quinolone" = "Quinolone", 
                         "rifamycin" = "Rifamycin",
                         "streptothricin" = "Streptothricin", 
                         "sulfonamide" = "Sulfonamide",
                         "tetracenomycin" = "Tetracenomycin",
                         "tetracycline" = "Tetracycline", 
                         "trimethoprim" = "Trimethoprim", 
                         "vancomycin" = "Vancomycin"))
unique_families <- unique(Gene_data_summarized$Family)
print(unique_families)


ARG_data_Add_summarized <- ARG_data_Add_summarized %>%
  mutate(Family = recode(Family,
                         "Antibacterial" = "Antibacterial_fatty_acid", 
                         "Beta-lactams" = "Beta_lactam",
                         "Pleuromutilin" = "Pleuromutilin_tiamulin"))

unique_families <- unique(ARG_data_Add_summarized$Family)
print(unique_families)

#####################################
### Step 4: Delete Families from Gene table that are not in ARG table 
#####################################

Gene_data_summarized <- Gene_data_summarized %>%
  filter(!Family %in% c("Bicyclomycin", "Bleomycin", "Factumycin", "other peptidic antibiotic", "Tetracenomycin" ))

ARG_data_Add_summarized <- ARG_data_Add_summarized %>%
  filter(!Family %in% c("Tunicamycin" ))


#####################################
### Step 5: Combine Tables 
#####################################

ARG_data_combined <- merge(ARG_data_Add_summarized, Gene_data_summarized, by = "Family", all = TRUE)
write.csv(ARG_data_combined, "Master/Metatranscriptomics/Excel_lists/ARG/ARG_data_combined.csv", row.names = FALSE)
