# Taxa Identification
The following script was used to download the viral reference database RefSeq from NCBI for Kraken2.
```bash
#!/bin/bash

# Define path to database
DB="/home/anna/KrakenLibraries"

#Download taxonomy data
#echo "Download taxonomy data..."
#kraken2-build --download-taxonomy --db $DB

#Download and set upviral sequence database (refseq)
echo "Set up viral Kraken2 database..."
kraken2-build --download-library viral --threads 12 --db $DB
#Create database
echo "Create Kraken2 database..."
kraken2-build --build --threads 12 --db $DB
#Check if database was successfully set up
if [ -f "$DB/hash.k2d" ]; then
    echo "Viral database was successfully downloaded and set up!"
else
    echo "Error in creating database!"
fi
```
With the next script the paired-end trimmed reads were taxonomically classified using the Kraken2 tool against the reference database RefSeq.
```bash
#!/bin/bash

#Set the directoies
KRAKEN2_PATH="/usr/lib/kraken2/kraken2"
DATABASE="/home/anna/KrakenLibraries/library/viral"
INPUTDIR="/media/2tbhdd/Antonia/Trimmed"
OUTPUTDIR="/media/2tbhdd/Antonia/kraken2_output"
WORKDIR="/home/anna/WorkingDir/Antonia"
ACCESSIONS_LIST="/home/anna/WorkingDir/Antonia/sra_accessions.txt"
LOGFILE="/home/anna/WorkingDir/Antonia/kraken2_process.log"

#Logfunction with timestamp
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$LOGFILE"
}

#Create the output directory if it does not exist
mkdir -p "$OUTPUTDIR"

#Process each accession
for accession in $(cat $ACCESSIONS_LIST); do
    log_message "Process $accession..."
    
    #Path to trimmed files
    FORWARD="$INPUTDIR/${accession}_1_val_1.fq"
    REVERSE="$INPUTDIR/${accession}_2_val_2.fq"
    
    #Check if trimmed files exist
    if [[ ! -f "$FORWARD" || ! -f "$REVERSE" ]]; then
        log_message "Error: One or both files for $accession missing!"
        continue
    fi

    #Move files into WorkingDir
    log_message "Move files in WorkingDIr..."
    mv "$FORWARD" "$WORKDIR"
    mv "$REVERSE" "$WORKDIR"
    wait
    log_message "Files successfully moved."
    
    #Start Kraken2
    log_message "Start Kraken2 for $accession..."
    "$KRAKEN2_PATH" \
        --db "$DATABASE" \
        --threads 4 \
        --report "$OUTPUTDIR/${accession}_kraken2_report.txt" \
        --output "$OUTPUTDIR/${accession}_kraken2_output.txt" \
        --unclassified-out "$OUTPUTDIR/${accession}_unclassified_#.fq" \
        --classified-out "$OUTPUTDIR/${accession}_classified_#.fq" \
        --paired "$WORKDIR/$(basename $FORWARD)" "$WORKDIR/$(basename $REVERSE)" &
    wait

    #Check results
    if [[ $? -eq 0 ]]; then
        log_message "Kraken2 for $accession successfully completed."
    else
        log_message "Error: Kraken2 not completed for $accession."
    fi
    
    #Move files back into InputDir
    log_message "Move files back into InputDir..."
    mv "$WORKDIR/$(basename $FORWARD)" "$INPUTDIR"
    mv "$WORKDIR/$(basename $REVERSE)" "$INPUTDIR"
done

#Done
log_message "Kraken2 process finished."
echo "All accessions were processed."
```
All Kraken2 report files were combined into a single table using the following R script.
```r
#####################################
### Load Required Libraries
#####################################
library(dplyr)

#Create list with all Kraken2 report files in directory
report_files <- list.files("/media/2tbhdd/Antonia/kraken2_output", pattern = "_kraken2_report.txt$", full.names = TRUE)

#Create empty dataframe to save all data
all_data <- data.frame()

#go through files
for (file in report_files) {
  #read data
  data <- read.table(file, header = FALSE, sep = "\t", stringsAsFactors = FALSE)
  
  #name columns
  colnames(data) <- c("Percentage", "Total_Reads", "Assigned_Reads", "Rank", "TaxID", "Taxon")
  
  #add column with accession numbers (SRR)
  srr_id <- sub(".*/(SRR[0-9]+)_kraken2_report.txt", "\\1", file)
  data$SRR_ID <- srr_id
  
  #combine data
  all_data <- bind_rows(all_data, data)
}

#Save data as CSV
combined_kraken2_report <- "/media/2tbhdd/Antonia/kraken2_output/combined_kraken2_report.csv"
write.csv(all_data, combined_kraken2_report, row.names = FALSE)
```
