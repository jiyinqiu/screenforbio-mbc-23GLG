#!/bin/bash

set -e
set -u
set -o pipefail


label=${1}
SCRIPTS=${2}

# get list of sp
grep ">" MIDORI_${label}.amp_blast.noN.edit.fasta | sed 's/>//g' | awk 'BEGIN{FS="_"}{print $1 "_" $2}' | sort -u > ${label}.sp_list.txt
# for each sp in list, if has >1 seq make temp fasta and collapse haplotypes, else print 1 seq
touch ${label}.uniq.fa
splist=($(cut -f1 ${label}.sp_list.txt))
for sp in ${splist[@]}
do
  grep ${sp} MIDORI_${label}.amp_blast.noN.edit.fasta | sed 's/>//g' | seqtk subseq MIDORI_${label}.amp_blast.noN.edit.fasta - > ${sp}.fa
  count=($(grep -c ">" ${sp}.fa))
  if [ $count -gt 1 ]
  then
    mafft --thread 4 --retree 2 --reorder ${sp}.fa > ${sp}.mafft.fa
    perl ${SCRIPTS}/collapsetypes_v4.6.pl -infile=${sp}.mafft.fa -nrdiffs=0
    seqbuddy ${sp}.mafft.fa.unique_haplotypes --clean_seq > ${sp}.uniq.fa
         #DY seqbuddy has a python error that adds "FTP Error: got more than 8192 bytes" to the end of ${sp}.uniq.fa
         #DY add this line to remove any instances of "FTP Error: got more than 8192 bytes"
         sed -i '/FTP Error: got more than 8192 bytes/d' ${sp}.uniq.fa
         sed -i '/FTP Error: \[Errno 8] nodename nor servname provided, or not known/d'  ${sp}.uniq.fa
    cat ${label}.uniq.fa ${sp}.uniq.fa > tmp
  else
    cat ${label}.uniq.fa ${sp}.fa > tmp
  fi
  mv tmp ${label}.uniq.fa
  rm ${sp}*
done
mv ${label}.uniq.fa MIDORI_${label}.amp_blast.noN.uniq.fasta
