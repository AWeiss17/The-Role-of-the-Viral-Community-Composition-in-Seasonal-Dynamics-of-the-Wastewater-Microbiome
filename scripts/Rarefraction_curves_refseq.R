#####################################
### Load Required Libraries
#####################################
library(ggplot2)
library(vegan)
library(dplyr)
library(scales)
library(readr)

#####################################
### Step 1: Load tables and set rownames (Transpose data)
#####################################

data_refseq_wide_SRR <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/count_data_refseq_wide_SRR.csv")


#####################################
### Step 2: Set Row Names & Clean Data
#####################################

data_refseq_wide_SRR <- as.data.frame(data_refseq_wide_SRR)

rownames(data_refseq_wide_SRR) <- data_refseq_wide_SRR$...1

data_refseq_wide_SRR <- data_refseq_wide_SRR[, -1]



#####################################
### Remove samples (Initial outliers from)
#####################################
# These 3 sample were removed from Nils in his initial Outlier screening
# The samples he identified as outliers are: SRR9006565, SRR9006522, SRR9006571
# Corresponding SRR_IDs for my dataset are: SRR9006560, SRR9006511, SRR9006558

removed_refseq_samples <- c("SRR9006560", "SRR9006511", "SRR9006558")
data_refseq_wide_SRR <- data_refseq_wide_SRR[!rownames(data_refseq_wide_SRR) %in% removed_refseq_samples, ]



#####################################
### Step 3: Rarefaction Curves
#####################################

rarecurves_refseq <- rarecurve(as.matrix(data_refseq_wide_SRR), step = 50)


# Assign each line in  plot to its corresponding AccessionID
names(rarecurves_refseq) <- rownames(data_refseq_wide_SRR)

# transform rarefaction data into data frame --> important for plotting
rare_data_refseq <- do.call(rbind, lapply(names(rarecurves_refseq), function(sample) { #names() = gives out names of the elements of the rare_curves data frame/table/list; then each is counted as the variable "sample"
  data.frame(
    Sample = sample,
    Reads = attr(rarecurves_refseq[[sample]], "Subsample"), #attr = Attribute of rare_curves(sample) (subsample: name given from rarefraction function)
    Family = rarecurves_refseq[[sample]]        # number of Genera per sample are taken from the rarecurve function/list for each 50 Read step
  )
}))



#####################################
### Step 4: Start Plotting Rarefraction curves
#####################################
# Create 48 distinct colors
num_samples <- length(unique(rare_data_refseq$Sample))  # Should be 48
color_palette <- hue_pal()(num_samples)  # Generates distinct colors

# Find label positions (last x,y value for each sample)
label_positions <- rare_data_refseq %>%
  group_by(Sample) %>%
  filter(Reads == max(Reads))

# Plot for refseq virus data

# Find dynamic x-axis and y-axis limits based on actual data
x_max <- max(rare_data_refseq$Reads, na.rm = TRUE)
y_max <- max(rare_data_refseq$Family, na.rm = TRUE)

# Create the plot  
ggplot(rare_data_refseq, aes(x = Reads, y = Family, group = Sample, color = Sample)) +
  geom_line(size = 1) +   
  scale_color_manual(values = color_palette) +   
  scale_x_continuous(limits = c(0, x_max * 1.05)) +  
  scale_y_continuous(limits = c(0, y_max * 1.05)) +  
  
  # Add lowest sequencing depth mindepth = 486
  geom_vline(xintercept = 486, linetype = "dashed", color = "black") +
  labs(
    title = "",
    x = "Number of Reads",
    y = "Number of Families",
    color = "Sample"  
  ) +
  
  theme_classic() +
  theme(
    legend.position = "right",  
    legend.key.height = unit(0.4, "cm"),  
    legend.key.width = unit(0.4, "cm"),  
    legend.text = element_text(size = 8),  
    legend.key.size = unit(0.4, "cm"),  
    legend.title = element_text(size = 9, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    
    panel.grid.major = element_line(size = 0.3, linetype = "dashed", color = "gray80"),
    panel.grid.minor = element_line(size = 0.2, linetype = "dotted", color = "gray90")
  ) +

  guides(color = guide_legend(ncol = 2))



ggsave(
  filename = "Rarefaction_Curve_RefSeq.png",  
  plot = last_plot(),  
  device = "png",  
  width = 10,  
  height = 7,  
  units = "in",  
  dpi = 600  
)




#####################################
### Step 5: Check read depth distribution
#####################################

# Summarize sequencing depths
seq_depths <- rowSums(data_refseq_wide_SRR)  # Sum of reads per sample

# Plot the distribution of sequencing depths
ggplot(data.frame(Sample = names(seq_depths), Depth = seq_depths), aes(x = Depth)) +
  geom_histogram(binwidth = 100, fill = "blue", color = "black", alpha = 0.7) +
  labs(title = "Sequencing Depth Distribution", x = "Number of Reads", y = "Number of Samples") +
  theme_minimal()



#####################################
### Step 6: Check how many samples would get lost at certain threshold
#####################################

# How many samples have lower reads than 500, 750, 1000?
table(rowSums(data_refseq_wide_SRR) < 400)
table(rowSums(data_refseq_wide_SRR) < 500)
table(rowSums(data_refseq_wide_SRR) < 600)
table(rowSums(data_refseq_wide_SRR) < 750)   
table(rowSums(data_refseq_wide_SRR) < 1000)  




###############################################################################
########################## Try different cut offs #############################
###############################################################################
# because sample with smallest read number is SRR9006549 with 486 reads (lowest sequencing depth)
# Sum all reads per sample
total_reads_per_sample <- rowSums(data_refseq_wide_SRR)

# find sample with lowest read count (sequencing depth)
min_reads_sample <- names(total_reads_per_sample)[which.min(total_reads_per_sample)]
min_reads_value <- min(total_reads_per_sample)

# Show sample with lowest read count (sequencing depth)
cat("Die Probe mit den wenigsten Reads ist:", min_reads_sample, "mit", min_reads_value, "Reads.\n")


#####################################
### Rarefy Data to Smallest Sequencing Depth
#####################################
# define sample with lowest read count (sequencing depth)
rarefy_depth <- min(rowSums(data_refseq_wide_SRR))

# Compute relative abundance per sample
data_relative <- data_refseq_wide_SRR / rowSums(data_refseq_wide_SRR)

# Scale to rarefied depth
data_rarefied <- round(data_relative * rarefy_depth)

# Remove taxa (columns) that are now all zeros
data_rarefied <- data_rarefied[, colSums(data_rarefied) > 0]

# Verify rarefaction result
summary(rowSums(data_rarefied))

# alternative approach for rarefying (with functionrrarefy)
#data_rarefied <- rrarefy(data_refseq_wide_SRR, sample = rarefy_depth)


#####################################
### Check if Taxa got Lost
#####################################

num_taxa_before <- ncol(data_refseq_wide_SRR)
num_taxa_after <- ncol(data_rarefied)

print(paste("Anzahl Taxa vorher:", num_taxa_before))
print(paste("Anzahl Taxa nach Rarefaction:", num_taxa_after))
print(paste("Verlorene Taxa:", num_taxa_before - num_taxa_after))
# no taxa got lost!

# Save the final rarefied table
write.csv(data_rarefied, paste0("Rarefied_count_data_refseq_486.csv"), row.names = TRUE, quote = FALSE)

