#####################################
### Load Required Libraries
#####################################
library(dplyr)

# Erstelle eine Liste aller Kraken2 Report-Dateien im Verzeichnis
report_files <- list.files("/media/2tbhdd/Antonia/kraken2_output", pattern = "_kraken2_report.txt$", full.names = TRUE)

# Leere DataFrame erstellen, um alle Daten zu speichern
all_data <- data.frame()

# Durchlaufe alle Dateien und lese sie ein
for (file in report_files) {
  # Lese die Datei ein
  data <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  # Benenne die Spalten
  colnames(data) <- c("Percentage", "Total_Reads", "Assigned_Reads", "Rank", "TaxID", "Taxon")
  
  # Füge eine neue Spalte hinzu, die den Dateinamen (SRR) enthält
  srr_id <- sub(".*/(SRR[0-9]+)_kraken2_report.txt", "\\1", file)
  data$SRR_ID <- srr_id
  
  # Füge die Daten dem Gesamtdatensatz hinzu
  all_data <- bind_rows(all_data, data)
}

# Speichere alle Daten als CSV-Datei
combined_kraken2_report <- "/media/2tbhdd/Antonia/kraken2_output/combined_kraken2_report.csv"
write.csv(all_data, combined_kraken2_report, row.names = FALSE)
