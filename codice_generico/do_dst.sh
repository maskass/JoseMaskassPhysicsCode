#! /bin/sh

ARG=1 # Numero di argomenti passati allo script

if [ $# -ne "$ARG" ]
# Verifica il numero degli argomenti e in caso ritorna errore
then
  echo "Inserire come argomento il nome della lista"
fi

file=$1
outdir=out_dst

mkdir -p $outdir

while read myline
do
outfile1="$(echo $myline | cut -d/ -f2)"
outfile2="$(echo $outfile1 | cut -d_ -f2)"

root -l -q macroCreateDst.cc+'("'$outfile2'")'
done < $file

exit 0
