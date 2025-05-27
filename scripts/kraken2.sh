#!/bin/bash

# Pfade zu den Verzeichnissen und Dateien
KRAKEN2_PATH="/usr/lib/kraken2/kraken2"
DATABASE="/home/anna/KrakenLibraries/library/viral"
INPUTDIR="/media/2tbhdd/Antonia/Trimmed"
OUTPUTDIR="/media/2tbhdd/Antonia/kraken2_output"
WORKDIR="/home/anna/WorkingDir/Antonia"
ACCESSIONS_LIST="/home/anna/WorkingDir/Antonia/sra_accessions.txt"
LOGFILE="/home/anna/WorkingDir/Antonia/kraken2_process.log"

# Logfunktion mit Timestamp
log_message() {
    echo "$(date '+%Y-%m-%d %H:%M:%S') - $1" >> "$LOGFILE"
}

# Erstelle das Ausgabeverzeichnis, falls es nicht existiert
mkdir -p "$OUTPUTDIR"

# Verarbeite jede Accession
for accession in $(cat $ACCESSIONS_LIST); do
    log_message "Verarbeite $accession..."
    
    # Pfade zu den getrimmten Dateien
    FORWARD="$INPUTDIR/${accession}_1_val_1.fq"
    REVERSE="$INPUTDIR/${accession}_2_val_2.fq"
    
    # Überprüfen, ob die getrimmten Dateien existieren
    if [[ ! -f "$FORWARD" || ! -f "$REVERSE" ]]; then
        log_message "FEHLER: Eine oder beide Dateien für $accession fehlen!"
        continue
    fi

    # Verschiebe Dateien ins Arbeitsverzeichnis
    log_message "Verschiebe die Dateien ins Arbeitsverzeichnis..."
    mv "$FORWARD" "$WORKDIR"
    mv "$REVERSE" "$WORKDIR"
    wait
    log_message "Dateien erfolgreich verschoben."
    
    # Starte Kraken2
    log_message "Starte Kraken2 für $accession..."
    "$KRAKEN2_PATH" \
        --db "$DATABASE" \
        --threads 4 \
        --report "$OUTPUTDIR/${accession}_kraken2_report.txt" \
        --output "$OUTPUTDIR/${accession}_kraken2_output.txt" \
        --unclassified-out "$OUTPUTDIR/${accession}_unclassified_#.fq" \
        --classified-out "$OUTPUTDIR/${accession}_classified_#.fq" \
        --paired "$WORKDIR/$(basename $FORWARD)" "$WORKDIR/$(basename $REVERSE)" &
    wait

    # Überprüfung der Ergebnisse
    if [[ $? -eq 0 ]]; then
        log_message "Kraken2 für $accession erfolgreich abgeschlossen."
    else
        log_message "FEHLER: Kraken2 konnte für $accession nicht abgeschlossen werden."
    fi
    
    # Verschiebe die Dateien zurück ins Ursprungsverzeichnis
    log_message "Verschiebe die Dateien zurück ins Ursprungsverzeichnis..."
    mv "$WORKDIR/$(basename $FORWARD)" "$INPUTDIR"
    mv "$WORKDIR/$(basename $REVERSE)" "$INPUTDIR"
done

# Abschlussnachricht
log_message "Kraken2-Prozess abgeschlossen."
echo "Alle Accessions wurden verarbeitet."
