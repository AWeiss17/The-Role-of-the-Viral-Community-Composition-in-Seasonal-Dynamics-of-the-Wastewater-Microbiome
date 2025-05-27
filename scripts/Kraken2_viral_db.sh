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
