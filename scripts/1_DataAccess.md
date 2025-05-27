# Data Access
Metagenomic read data in FASTQ format were obtained from the Sequence Read Archive (SRA). The data download was performed using the prefetch command to retrieve the corresponding SRA files. The following script was used to download the SRA files.
## Prefetch

The following script was used to download FASTQ files from the SRA using the `prefetch` command. The list oa accessions `(sra_accession.txt)` can be found in the `Data` folder of this repository.

```bash
#!/bin/bash
#Prefetches from SRA by accession stored in textfile "sra_accessions.txt" and by number of accessions in that file.
#Takes two arguments: the line number of the first and the last accession in #that file to prefetch [1-48].
#Go through all accessions.
#"echo" to get information about the current process
i=1
while read -r accession; do
  echo "looking at $accession"
  # the prefetch command
  /home/anna/sratoolkit.3.0.7-ubuntu64/bin/prefetch "$accession" -O /home/anna/WorkingDir/Antonia &
  wait
  echo "prefetched $accession"
  # move the sra file to the 2tbhdd drive
  mv "/home/anna/WorkingDir/Antonia/$accession/$accession.sra" "/media/2tbhdd/Antonia" &
  wait
  echo "moved $accession"
  # look at the next accession
  let "i=i+1"
done < sra_accessions.txt
```

With the following sra.getter.sh script the SRA files were moved from the external hard drive to the local working directory and back to their original location.
```bash
#!/bin/bash
#move the sra file called by its accession in the first argument back where it can be used (it was stored on the 2TB drive, that is echoed.)
#second argument is either "get" or "putback". "Get" outputs a string, which can be received as third argument in "putback" mode to put it back to its original saving position.
accession="$1"
operation="$2"
knownlocation="$3"
if [ "$operation" = "get" ]
then
  if [ $(ls /media/2tbhdd/Antonia/ | grep "$accession"\.sra | wc -l) -gt 0 ]
  then
   mv /media/2tbhdd/Antonia/"$accession".sra /home/anna/WorkingDir/Antonia/$accession &
   wait
   echo tbd
  else
   echo wdd
  fi
fi
if [ "$operation" = "putback" ] && [ "$knownlocation" = "tbd" ]
then
mv "/home/anna/WorkingDir/Antonia/$accession/$accession.sra" "/media/2tbhdd/Antonia" &
wait
fi
```
