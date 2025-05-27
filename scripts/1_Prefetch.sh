#!/bin/bash
#prefetches from SRA by accession stored in textfile "sra_accessions.txt" and by number of accessions in that
#file.
#takes two arguments: the line number of the first and the last accession in that file to prefetch.
# go through all accession
#"echo" is just to get information about the currend process
i=1
while read -r accession; do
echo "looking at $accession"
# the prefetch command
/home/anna/sratoolkit.3.0.7-ubuntu64/bin/prefetch "$accession" -O /home/anna/WorkingDir/Antonia &
wait
echo "prefetched $accession"
# move the sra file to the terabyte drive
mv "/home/anna/WorkingDir/Antonia/$accession/$accession.sra" "/media/2tbhdd/Antonia" &
wait
echo "moved $accession"
# look at the next accession
let "i=i+1"
done < sra_accessions.txt
