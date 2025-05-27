# Statistical Analyses
- `ViolinPlots.R`: Violin plots were generated to provide a first overview of the distribution of raw read counts per sample across seasons.
- `RelativeAbundance.R`: Relative abundance plots were generated for viral, prokaryotic, eukaryotic and ARG data to visualize the temporal distribution of the most abundant taxa across all samples.
- `NMDS_Viruses/Prokaryotes/Eukaryotes/ARGs.R`: Non-metric multidimensional scaling (NMDS) was performed with rarefied countdata to visualize overall community structure and seasonal patterns. To explore potential associations, environmental variables, as well as the most abundant taxa of viruses, prokaryotes, eukaryotes, and ARGs were fitted as vectors onto the ordination space using the envfit() function.
- `PCoA.R`: with this script Principal Coordinate Analysis (PCoA) based on Bray-Curtis dissimilarities were performed based on rarefied viral, prokaryotic, eukaryotic, and ARG count data.
- `PERMANOVA.R`: PERMANOVA analyses were performed to assess if environmental variables, seasonal variation, and the community composition of viruses, prokaryotes, eukaryotes, and ARGs have a significant influence on the beta diversity of viral, prokaryotic and ARG communities. For these calculations, PCoA axes explaining the highest variance were included as explanatory variables.
- `SpearmanCorrelation.R`:
- `SEM.R`: SEM was applied using Shannon diversity indices to assess seasonal effects and associations between viral, prokaryotic, eukaryotic, and ARG communities. The model was specified in lavaan and evaluated based on unstandardized path coefficients.

