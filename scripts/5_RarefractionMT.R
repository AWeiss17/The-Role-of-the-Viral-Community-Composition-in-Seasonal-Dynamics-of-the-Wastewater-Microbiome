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

data_pr2_wide_SRR <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/count_data_pr2_wide_SRR.csv")
data_silva_wide_SRR <- read_csv("Master/Metatranscriptomics/Excel_lists/Meta_omic Tables/count_data_silva_wide_SRR.csv")

#####################################
### Step 2: Combine PR2 and SILVA Data
#####################################

data_RNA <- merge(
  data_pr2_wide_SRR, data_silva_wide_SRR, 
  by = "...1",  # Merge by first column (sample ID)
  all = TRUE    # Keep all rows, fill missing values with NA
)

#####################################
### Step 3: Set Row Names & Clean Data
#####################################

data_RNA <- as.data.frame(data_RNA)

rownames(data_RNA) <- data_RNA$...1

data_RNA <- data_RNA[, -1]  # Remove first column after setting row names


#####################################
### Remove samples (Initial outliers from)
#####################################
# These 3 sample were removed from Nils in his initial Outlier screening
# The samples he identified as outliers are: SRR9006565, SRR9006522, SRR9006571
# Corresponding SRR_IDs for my dataset are: SRR9006560, SRR9006511, SRR9006558

removed_RNA_samples <- c("SRR9006565", "SRR9006522", "SRR9006571")
data_RNA <- data_RNA[!rownames(data_RNA) %in% removed_RNA_samples, ]


#####################################
### Step 4: Rarefaction Curves
#####################################

rarecurves_RNA <- rarecurve(as.matrix(data_RNA), step = 50)


# Assign each line in  plot to its corresponding AccessionID
names(rarecurves_RNA) <- rownames(data_RNA)

# transform rarefaction data into data frame --> important for plotting
rare_data_RNA <- do.call(rbind, lapply(names(rarecurves_RNA), function(sample) { #names() = gives out names of the elements of the rare_curves data frame/table/list; then each is counted as the variable "sample"
  data.frame(
    Sample = sample,
    Reads = attr(rarecurves_RNA[[sample]], "Subsample"), #attr = Attribute of rare_curves(sample) (subsample: name given from rarefraction function)
    Genera = rarecurves_RNA[[sample]]        # number of Genera per sample are taken from the rarecurve function/list for each 50 Read step
  )
}))



#####################################
### Step 5: Start Plotting Rarefraction curves
#####################################
# Create 48 distinct colors
num_samples <- length(unique(rare_data_RNA$Sample))  # Should be 48
color_palette <- hue_pal()(num_samples)  # Generates distinct colors

# Find label positions (last x,y value for each sample)
label_positions <- rare_data_RNA %>%
  group_by(Sample) %>%
  filter(Reads == max(Reads))

# Plot for Metatranscriptomic RNA data (Nils)

# Find dynamic x-axis and y-axis limits based on actual data
x_max <- max(rare_data_RNA$Reads, na.rm = TRUE)
y_max <- max(rare_data_RNA$Genera, na.rm = TRUE)

# Create the plot  
ggplot(rare_data_RNA, aes(x = Reads, y = Genera, group = Sample, color = Sample)) +
  geom_line(size = 1) +   # Increase line thickness
  scale_color_manual(values = color_palette) +   # Apply distinct colors
  # Add cutt off at mindepth = 1e6
  geom_vline(xintercept = 1e6, linetype = "dashed", color = "black") +
  # Labels & Titles
  labs(
    title = "Rarefaction Curves",
    x = "Number of Reads",
    y = "Number of Genera",
    color = "Sample Names"  # Legend title
  ) +
  
  # Improve Theme  
  theme_classic() +
  theme(
    legend.position = "right",  # Keep legend for color reference
    legend.key.height = unit(0.4, "cm"),  # Reduce spacing between legend items
    legend.key.width = unit(0.4, "cm"),  # Reduce width of legend keys
    legend.text = element_text(size = 8),  # Adjust legend text size
    legend.key.size = unit(0.4, "cm"),  # Reduce legend spacing
    legend.title = element_text(size = 9, face = "bold"),
    plot.title = element_text(size = 14, face = "bold", hjust = 0.5),  # Center title
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12, face = "bold"),
    
    # Add Grid Lines
    panel.grid.major = element_line(size = 0.3, linetype = "dashed", color = "gray80"),
    panel.grid.minor = element_line(size = 0.2, linetype = "dotted", color = "gray90")
  ) +
  # Make legend two columns
  guides(color = guide_legend(ncol = 2))



ggsave(
  filename = "Rarefaction_Curve_RNA_1e6_48samples.png",  # Output file name
  plot = last_plot(),  # Save the last generated plot
  device = "png",  # Save as PNG
  width = 10,  # Adjust width (in inches)
  height = 7,  # Adjust height (in inches)
  units = "in",  # Units for width and height
  dpi = 600  # High resolution (300-600 DPI recommended for publication)
)


#####################################
### Step 6: Cut off at min_depth = 1e6
#####################################
# at a cut off at 1000000 reads
mindepth <- 1000000

# Check which sample get lost
removed_samples <- rownames(data_RNA)[rowSums(data_RNA) <= mindepth]
print(paste0("Removed sample(s): ", paste(removed_samples, collapse=", ")))  

# Filter out samples below cutoff
data_filtered <- data_RNA[rowSums(data_RNA) >= mindepth, ]

# Show distribution of Seq. depth
summary(rowSums(data_filtered))

hist(rowSums(data_filtered), breaks = 20, col = "lightblue", main = "Distribution of sequencing depths", xlab = "Reads per sample")


# Find the smallest remaining sequencing depth (target for rarefaction)
rarefy_depth <- min(rowSums(data_filtered))  

print(paste("Rarefaction target depth:", rarefy_depth))  



#####################################
### Step 7: Rarefy Data to Smallest Remaining Depth = 1066681
#####################################

# Compute relative abundance per sample
data_relative <- data_filtered / rowSums(data_filtered)

# Scale to rarefied depth
data_rarefied <- round(data_relative * rarefy_depth)

# Remove taxa (columns) that are now all zeros
data_rarefied <- data_rarefied[, colSums(data_rarefied) > 0]

# Verify rarefaction result
summary(rowSums(data_rarefied)) 


# Save the final rarefied table
write.csv(data_rarefied, paste0("Rarefied_count_data_RNA.csv"), row.names = TRUE, quote = FALSE)

#####################################
### Step 9: Check if Taxa got Lost
#####################################

num_taxa_before <- ncol(data_filtered)
num_taxa_after <- ncol(data_rarefied)

print(paste("Anzahl Taxa vorher:", num_taxa_before))
print(paste("Anzahl Taxa nach Rarefaction:", num_taxa_after))
print(paste("Verlorene Taxa:", num_taxa_before - num_taxa_after))

# Show which taxa were removed
removed_taxa <- setdiff(colnames(data_filtered), colnames(data_rarefied))
print("Entfernte Taxa:")
print(removed_taxa)
      
