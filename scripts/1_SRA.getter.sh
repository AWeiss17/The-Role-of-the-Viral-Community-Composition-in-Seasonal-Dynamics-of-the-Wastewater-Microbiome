#!/bin/bash
# move the sra file called by its accession in the first argument back where it can be used (it was stored on the
#2TB drive, that is echoed.)
# second argument is either "get" or "putback". "Get" outputs a string, which can be received as third argument in
#"putback" mode to put it back to its original saving position.
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
