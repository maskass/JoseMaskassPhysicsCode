#! /bin/sh

ARG=1 # Numero di argomenti passati allo script

E_ERR_ARG=65   # codice di errore

if [ $# -ne "$ARG" ]
# Verifica il numero degli argomenti e in caso ritorna errore
then
  echo "Inserire come argomento il nome della lista"
  exit $E_ERR_ARG
fi

file=$1

mkdir -p root

while read myline
do
outfile1="$(echo $myline | cut -d/ -f2)"
outfile2="$(echo $outfile1 | cut -d. -f1)"
outfile3="$(echo $outfile2 | cut -d_ -f1)"
outfile4="$(echo $outfile2 | cut -d_ -f3)"

root -l -q macro.cc'("'$outfile3'","'$outfile4'")'

ls root/$outfile3* > $outfile3.list

done < $file

exit 0
