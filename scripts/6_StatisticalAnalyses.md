# Statistical Analyses
- `ViolinPlots.R`: Violin plots were generated to provide a first overview of the distribution of raw read counts per sample across seasons.
- `RelativeAbundance.R`: Relative abundance plots were generated for viral, prokaryotic, eukaryotic and ARG data to visualize the temporal distribution of the most abundant taxa across all samples.
- `NMDS_Viruses/Prokaryotes/Eukaryotes/ARGs.R`: Non-metric multidimensional scaling (NMDS) was performed with rarefied countdata to visualize overall community structure and seasonal patterns. To explore potential associations, environmental variables, as well as the most abundant taxa of viruses, prokaryotes, eukaryotes, and ARGs were fitted as vectors onto the ordination space using the envfit() function.
