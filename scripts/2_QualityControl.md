# Quality Control
With the command fasterq-dump, the SRA files were converted into FASTQ format, followed by quality control with FastQC.
```bash
#!/bin/bash
#make FASTQs from SRA by accession stored in "sra_accessions.txt" textfile
#takes two arguments: the line number of the first accession in #"sra_accessions.txt" textfile
#and the last accession in that file to make FASTQs from [1-48]
#perform FASTQC quality control and delete the FASTQs
#go through all accessions
i=1
for accession in $(cat sra_accessions.txt)
do
 # look only at selected accessions
if [ "$i" -ge "$1" ] && [ "$i" -le "$2" ]
then
  echo "looking at $accession"
  bash sra.getter.sh "$accession" get > "${accession}_tempfile.txt"
  /home/anna/sratoolkit.3.0.7-ubuntu64/bin/fasterq-dump -S "$accession"
  bash sra.getter.sh "$accession" putback "$(cat "${accession}_tempfile.txt")"
  rm "${accession}_tempfile.txt"
# Check if the fastq files are in the external drive, and move them back if needed
 if [ -f "/media/2tbhdd/Antonia/${accession}_1.fastq" ]; then
   mv "/media/2tbhdd/Antonia/${accession}_1.fastq" "/home/anna/WorkingDir/Antonia/"
 fi
 if [ -f "/media/2tbhdd/Antonia/${accession}_2.fastq" ]; then
   mv "/media/2tbhdd/Antonia/${accession}_2.fastq" "/home/anna/WorkingDir/Antonia/"
 fi
   fastqc "/home/anna/WorkingDir/Antonia/${accession}_1.fastq" "/home/anna/WorkingDir/Antonia/${accession}_2.fastq" -q &
   wait
 if [ -f "/home/anna/WorkingDir/Antonia/${accession}_1.fastq" ]; then
   mv "/home/anna/WorkingDir/Antonia/${accession}_1.fastq" "/media/2tbhdd/Antonia"
 fi
 if [ -f "/home/anna/WorkingDir/Antonia/${accession}_2.fastq" ]; then
   mv "/home/anna/WorkingDir/Antonia/${accession}_2.fastq" "/media/2tbhdd/Antonia"
 fi
   wait
#rm "${accession}"*.fastq
fi
# look at the next accession
((i=i+1))
done
```
