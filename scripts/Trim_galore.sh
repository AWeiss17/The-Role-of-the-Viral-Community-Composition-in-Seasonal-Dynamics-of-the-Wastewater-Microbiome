#!/bin/bash
#Makes FASTQs from SRA by accession, then trimming of FASTQs
#Set the directories
OUTPUTDIR_TRIM="/media/2tbhdd/Antonia/Trimmed"
INPUTDIR="/media/2tbhdd/Antonia"
WORKDIR="/home/anna/WorkingDir/Antonia"
ACCESSIONS_LIST="/home/anna/WorkingDir/Antonia/sra_accessions.txt"
LOGFILE="/home/anna/WorkingDir/Antonia/trim_process.log"

# Log function with timestamp
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" >> $LOGFILE
}

#Go through all accessions in sra_accessions.txt
i=1
for accession in $(cat $ACCESSIONS_LIST); do
    log_message "Process $accession..."

    #move original FASTQ files
    FORWARD_FASTQ="${INPUTDIR}/${accession}_1.fastq"
    REVERSE_FASTQ="${INPUTDIR}/${accession}_2.fastq"

    if [ -f "$FORWARD_FASTQ" ] && [ -f "$REVERSE_FASTQ" ]; then
        log_message "Move $FORWARD_FASTQ and $REVERSE_FASTQ in WorkingDir..."
        mv "$FORWARD_FASTQ" "$WORKDIR" &
        mv "$REVERSE_FASTQ" "$WORKDIR" &
        wait
        log_message "Files successfully moved."
    else
        log_message "Error: one or both FASTQ files for $accession missing"
    fi

    #process Trim Galore 
    if [ -f "$WORKDIR/${accession}_1.fastq" ] && [ -f "$WORKDIR/${accession}_2.fastq" ]; then
        log_message "Start Trimming-Process for $accession..."
        /usr/bin/trim_galore --paired -q 30 --three_prime_clip_R1 5 --three_prime_clip_R2 5 --output_dir "$WORKDIR" "$WORKDIR/${accession}_1.fastq" "$WORKDIR/${accession}_2.fastq" &
        wait
        log_message "Trimming for $accession completed."
    else
        log_message "Error: one or both FASTQ files for trimming are missing!"
    fi

    #move not trimmed FASTQ files back
    if [ -f "$WORKDIR/${accession}_1.fastq" ] && [ -f "$WORKDIR/${accession}_2.fastq" ]; then
        log_message "Move not trimmed FASTQ files back in OUTPUTDIR..."
        mv "$WORKDIR/${accession}_1.fastq" "$INPUTDIR" &
        mv "$WORKDIR/${accession}_2.fastq" "$INPUTDIR" &
        wait
        log_message "Not trimmed files successfully moved."
    else
        log_message "Erorr: not trimmed files missing!"
    fi

    #move trimmed files
    if [ -f "$WORKDIR/${accession}_1_val_1.fq" ] && [ -f "$WORKDIR/${accession}_2_val_2.fq" ]; then
        log_message "Move trimmed files back in OUTPUTDIR_TRIM..."
        mv "$WORKDIR/${accession}_1_val_1.fq" "$OUTPUTDIR_TRIM" &
        mv "$WORKDIR/${accession}_2_val_2.fq" "$OUTPUTDIR_TRIM" &
        wait
        log_message "Trimmed files successfully moved."
    else
        log_message "Error: trimmed files missing!"
    fi

    log_message "$accession processed!"
done

# Done
log_message "Trimming-process done. All files were moved back."
