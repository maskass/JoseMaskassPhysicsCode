#! /bin/sh

ARG=1 # Numero di argomenti passati allo script

if [ $# -ne "$ARG" ]
# Verifica il numero degli argomenti e in caso ritorna errore
then
  echo "Inserire come argomento il nome della lista"
fi

file=$1

while read myline
do
mimerge $myline
done < $file

exit 0
